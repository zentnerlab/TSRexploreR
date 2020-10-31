
#' Import BAMs
#'
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#'
#' @param experiment tsr explorer object
#' @param paired Whether the BAMs are paired (TRUE) or unpaired (FALSE)
#' @param sample sheet to import bams
#' @param soft_remove Remove read if greater than this number of soft-clipped bases.
#' @param proper_pair Whether reads should be properly paired for paired end data.
#'   TRUE by default when data is paired end.
#' @param remove_seconday Remove non-primary reads.
#'
#' @export

import_bams <- function(
  experiment,
  paired,
  sample_sheet=NULL,
  soft_remove=3,
  proper_pair=NULL,
  remove_secondary=TRUE
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is.null(sample_sheet) ||
    (is.character(sample_sheet) | is.data.frame(sample_sheet))
  )
  assert_that(is.null(soft_remove) || is.count(soft_remove))
  assert_that(is.null(proper_pair) || is.flag(proper_pair))
  assert_that(is.flag(remove_secondary))

  ## Prepare sample sheet if required.
  sample_sheet_type <- case_when(
    is.null(sample_sheet) ~ "internal",
    is.character(sample_sheet) ~ "file",
    is.data.frame(sample_sheet) ~ "data_frame"
  )

  sample_sheet <- switch(sample_sheet_type,
    "internal"=experiment@meta_data$sample_sheet,
    "file"=fread(sample_sheet, sep="\t"),
    "data_frame"=as.data.table(sample_sheet)
  )

  samples <- as.list(sample_sheet[, file_1])
  names(samples) <- sample_sheet[, sample_name]

  ## Specifying flag settings for filtering.
  flag_args <- list(
    isUnmappedQuery=FALSE # Remove unmapped.
  )
  if (remove_secondary) {
    flag_args <- c(flag_args, list(
      isSecondaryAlignment=FALSE # Remove non-primary.
    ))
  }
  if (paired & (!is.null(proper_pair) && proper_pair)) {
    flag_args <- c(flag_args, list(
      isPaired=TRUE, # Return only paired reads.
      isProperPair=TRUE, #Remove improperly paired reads.
      hasUnmappedMate=FALSE #Remove reads with unmapped mate.
    ))
  }

  ## Import BAMs.
  if (paired) {
    bams <- map(samples, function(x) {
      bam <- readGAlignmentPairs(x, param=ScanBamParam(what="seq", flag=do.call(scanBamFlag, flag_args)))
      bam <- as.data.table(bam)[, .(
        seqnames=seqnames.first, start=start.first, end=end.first,
        strand=strand.first, cigar=cigar.first, seq=seq.first
      )]
      return(bam)
    })
  } else {
    bams <- map(samples, function(x) {
      bam <- readGAlignments(x, param=ScanBamParam(what="seq", flag=do.call(scanBamFlag, flag_args)))
      bam <- as.data.table(bam)[, .(seqnames, start, end, strand, cigar, seq)]
      return(bam)
    })
  }

  ## Get soft-clipped bases.
  walk(bams, function(x) {
    x[, n_soft := as.numeric(ifelse(
      strand == "+",
      str_extract(cigar, "^[[:digit:]]+(?=S)"), # For + strand soft-clip is at cigar beginning.
      str_extract(cigar, "[[:digit:]]+S$")      # For - strand soft-clip is at cigar end.
    ))][,
      n_soft := replace_na(n_soft, 0)
    ][
      n_soft == 0, seq := "none" # If there are no soft-clipped bases, replace seq with 'none'
    ][
      n_soft > 0, # Only subset sequences with soft-clipped bases.
      seq := str_sub(seq, end=n_soft) # The sequence for negative strand is reverse complement.
    ]
    setnames(x, old="seq", new="seq_soft")
    return(x)
  })

  ## Remove reads with too many soft-clipped bases.
  bams <- map(bams, ~ .x[is.na(n_soft) | n_soft <= soft_remove])

  ## Convert to GRanges.
  walk(bams, function(x) {
    x[,
      start := ifelse(strand == "+", start, end)
    ][,
      end := start
    ]
    x[, cigar := NULL]
  })

  bams <- map(bams, as_granges)

  ## Add GRanges to tsr explorer object.
  experiment@experiment$TSSs <- bams

  return(experiment)
}

#' Aggregate TSSs
#'
#' @description
#' Aggregate overlapping TSSs into a total sum score.
#'
#' @param experiment tsr explorer object
#'
#' @export

tss_aggregate <- function(experiment) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))

  ## Pull samples out.
  samples <- experiment@experiment$TSSs

  ## Aggregate TSSs.
  samples <- map(samples, function(x) {
    x <- as.data.table(x)
    x <- x[, .(score=.N), by=.(seqnames, start, end, strand)]
    x <- as_granges(x)
    return(x)
  })

  ## Add samples back.
  experiment@experiment$TSSs <- samples

  return(experiment)
}

#' Correct G
#'
#' @description
#' Correct overrepresentation of G TSSs
#'
#' @param experiment tsr explorer object
#' @param assembly Fasta assembly of genome or BSgenome object.
#'
#' @export

G_correction <- function(
  experiment,
  assembly
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(assembly) | is(assembly, "BSgenome"))

  ## Prepare assembly.
  assembly_type <- case_when(
    is.character(assembly) ~ "file",
    is(assembly, "BSgenome") ~ "bsgenome"
  )

  assembly <- switch(
    assembly_type,
    "file"=FaFile(assembly),
    "bsgenome"=assembly
  )

  ## Get samples.
  select_samples <- experiment@experiment$TSSs

  ## Retrieve +1 base.
  select_samples <- map(select_samples, function(x) {
    x$plus_one <- as.character(getSeq(assembly, x))
    x <- as.data.table(x)
    return(x)
  })

  ## Correct for frequency of soft-clipped Gs and templated Gs.
  select_samples <- map(select_samples, function(x) {

    # Frequency of soft-clipped Gs.
    sfreq <- x[, .(count=.N), by=seq_soft]
    sfreq[, freq := count / sum(count)]
    sfreq <- sfreq[seq_soft == "G", freq]

    # Correct frequencies.
    x[
      plus_one == "G" & n_soft == 0,
      start := ifelse(
        rbinom(nrow(x[plus_one == "G" & n_soft == 0]), 1, sfreq) == 1,
        ifelse(strand == "+", start + 1, start - 1),
        start
      )
    ][,
      c("end", "plus_one") := list(start, NULL)
    ]

    return(x)
  })

  ## Add data back to tsr explorer object.
  select_samples <- map(select_samples, as_grages)
  experiment@experiment$TSSs <- select_samples

  return(experiment)
}
