#' Import BAMs
#'
#' @description 
#' Import BAM files with optional quality control parameters.
#'
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments
#'
#' @inheritParams common_params
#' @param paired Whether the BAMs are paired (TRUE) or unpaired (FALSE).
#' @param soft_remove Remove read if greater than this number of soft-clipped bases 
#'   are present at its 5' most end.
#' @param proper_pair Remove reads without a proper pair SAM flag.
#'   TRUE by default when data is paired-end.
#' @param remove_seconday Remove reads with non-primary SAM flag set (TRUE).
#' @param remove_duplicate Remove reads with duplicate SAM flag set (TRUE).
#'
#' @details
#' Import BAMs using the information from the sample sheet.
#' If the BAMs are from paired end data,
#'   'proper_pair' allows removal of reads without a proper-pair SAM flag.
#' Additionally 'remove_seconday' and 'remove_duplicate' will remove reads
#'   with the secondary alignments and duplicate flags set.
#'
#' Most TSS mapping methodologies tend to add at least one non-templated base
#'   at the 5' end of the read.
#' Futhermore, template switching reverse transcription (TSRT) methods such
#'   as STRIPE-seq or nanoCAGE can have up to 3 or 4 non-templated 5' bases.
#' We recommend setting `soft_remove` to at minimum 3 because of this,
#'   Which removes the read if the given number of soft-clip bases are exceeded.
#'
#' @return TSRexploreR object with BAM GRanges and soft-clip information.
#'
#' @examples
#' bam_file <- system.file("extdata", "S288C.bam", package="TSRexploreR")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#' samples <- data.frame(sample_name="S288C", file_1=bam_file, file_2=NA)
#'
#' bam_file %>%
#'   tsr_explorer(sample_sheet=samples, genome_assembly=assembly) %>%
#'   import_bams(paired=TRUE)
#'
#' @export

import_bams <- function(
  experiment,
  paired,
  sample_sheet=NULL,
  soft_remove=3,
  proper_pair=NULL,
  remove_secondary=TRUE,
  remove_duplicate=FALSE
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.flag(paired))
  assert_that(
    is.null(sample_sheet) ||
    (is.character(sample_sheet) | is.data.frame(sample_sheet))
  )
  assert_that(is.null(soft_remove) || is.count(soft_remove))
  assert_that(is.null(proper_pair) || is.flag(proper_pair))
  assert_that(is.flag(remove_secondary))
  assert_that(is.flag(remove_duplicate))

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
  if (remove_duplicate) {
    flag_args <- c(flag_args, list(isDuplicate=FALSE))
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

  ## Add GRanges and sample sheet to tsr explorer object.
  experiment@experiment$TSSs <- bams
  experiment@meta_data$sample_sheet <- sample_sheet

  return(experiment)
}

#' Aggregate TSSs
#'
#' @description
#' Aggregate overlapping TSSs into a total sum score.
#'
#' @inheritParams common_params
#'
#' @export

tss_aggregate <- function(experiment) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))

  ## Get samples.
  samples <- experiment@experiment$TSSs

  ## Aggregate TSSs.
  samples <- map(samples, function(x) {
    x <- as.data.table(x)
    x <- x[, .(score=.N), by=.(seqnames, start, end, strand)]
    x <- as_granges(x)
    return(x)
  })

  ## Add samples back to TSRexploreR object.
  experiment@experiment$TSSs <- samples

  return(experiment)
}

#' Correct G artifact
#'
#' @description
#' Correct overrepresentation of 5' G bases added during reverse transcription.
#'
#' @param experiment TSRexploreR object.
#' @param assembly Genome assembly in FASTA or BSgenome format.
#'
#' @export

G_correction <- function(
  experiment,
  assembly=NULL
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is.null(assembly) ||
    (is.character(assembly) | is(assembly, "BSgenome"))
  )

  ## Prepare assembly.
  assembly <- .prepare_assembly(assembly, experiment)

  assembly_type <- case_when(
    is(assembly, "FaFile") ~ "fafile",
    is(assembly, "BSgenome") ~ "bsgenome"
  )

  ## Get samples.
  select_samples <- experiment@experiment$TSSs

  ## Retrieve +1 base.
  select_samples <- map(select_samples, function(x) {
    seq <- switch(
      assembly_type,
      "fafile"=Rsamtools::getSeq(assembly, x),
      "bsgenome"=BSgenome::getSeq(assembly, x)
    )
    x$plus_one <- as.character(seq)
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

  ## Add data back to TSRexploreR object.
  select_samples <- map(select_samples, as_granges)
  experiment@experiment$TSSs <- select_samples

  return(experiment)
}
