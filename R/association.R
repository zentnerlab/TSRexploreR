#' Merge Samples
#'
#' Merge replicates or selected samples.
#'
#' @param experiment tsrexplorer object
#' @param data_type Whether to merge TSSs or TSRs
#' @param merge_replicates If 'TRUE', replicate groups will be merged
#' @param threshold Filter out TSSs or TSRs below this raw count threshold before merging
#' @param sample_list If merge_replicates is set to 'FALSE',
#' specify what samples to merge in list format.
#'
#' @rdname merge_samples-function
#' @export

merge_samples <- function(
  experiment,
  data_type = c("tss", "tsr"),
  threshold = NA,
  merge_replicates = FALSE,
  sample_list = NA
) {
  
  ## Prepare list of samples to be merged.
  if (merge_replicates) {
    data_column <- ifelse(data_type == "tss", "tss_name", "tsr_name")

    sample_list <- experiment@meta_data$sample_sheet[type == data_type] %>%
      split(.$replicate_id) %>%
      map(~ pull(., name))
  }

  ## Merge feature sets.
  if (data_type == "tss") {
    merged_samples <- map(sample_list, function(sample_group) {
      select_samples <- experiment@experiment$TSSs[sample_group]
      select_samples <- map(select_samples, as.data.table)

      merged <- rbindlist(select_samples)
      if (!is.na(threshold)) merged <- merged[score >= threshold]
      merged <- merged[, .(score = sum(score)), by = .(seqnames, start, end, strand)]

      return(merged)
    })
  } else if (data_type == "tsr") {
    merged_samples <- map(sample_list, function(sample_group) {

      select_samples <- experiment@experiment$TSRs[sample_group]
      if (!is.na(threshold)) {
        select_samples <- map(select_samples, function(x) {
          x <- x[score(x) >= threshold]
          return(x)
        })
      }

      tsr_consensus <- select_samples %>%
        purrr::reduce(c) %>%
        GenomicRanges::reduce(ignore.strand = FALSE)

      merged <- map(select_samples, function(x) {
          overlap <- findOverlapPairs(query = tsr_consensus, subject = x)
          overlap <- as.data.table(overlap)
      })
      merged <- rbindlist(merged)[,
        .(first.seqnames, first.start, first.end,
        first.strand, second.X.score)
      ]

      setnames(
        merged,
        old = c(
          "first.seqnames", "first.start", "first.end",
          "first.strand", "second.X.score"
        ),
        new = c("seqnames", "start", "end", "strand", "score")
      )

      merged <- merged[,
        .(score = sum(score)),
        by = .(seqnames, start, end, strand)
      ]

      return(merged)
    })
  }

  ## Convert merged samples to GRanges.
  merged_samples <- map(merged_samples, as_granges)

  ## Return merged samples.
  if (data_type == "tss") {
    experiment@experiment$TSSs <- c(experiment@experiment$TSSs, merged_samples)
  } else if (data_type == "tsr") {
    experiment@experiment$TSRs <- c(experiment@experiment$TSRs, merged_samples)
  }

  return(experiment)
}

#' Associate TSSs
#'
#' @description
#' Associate TSSs with TSRs
#'
#' @importFrom plyranges join_overlap_left_directed
#'
#' @param experiment tsrexplorer object
#' @param use_sample_sheet Whether to use a sample sheet as a key for association of TSS and TSR samples
#' @param sample_list If 'use_sample_sheet' is FALSE, provide a list with TSR and TSS sample names
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#' tsre_exp <- associate_with_tsr(tsre_exp, sample_list = list(
#'   "S288C_WT_1" = "S288C_WT_1", "S288C_WT_2" = "S288C_WT_2", "S288C_WT_3" = "S288C_WT_3",
#'   "S288C_D_1" = "S288C_D_1", "S288C_D_2" = "S288C_D_2", "S288C_D_3" = "S288C_D_3"
#' ))
#'
#' @details
#' TSRexploreR provides many options for the import and merging of TSSs and TSRs.
#' Because of this, TSS samples must be associated with TSR samples after TSR import,
#'   TSR merging, or TSS clustering using this function.
#' This adds an extra workflow step, but provides more analytical flexibility (qq meaning what?).
#' Each TSS with genomic coordinates overlapping those of a TSR in the specified TSR sample
#'   will be linked to that TSR.
#' TSSs not overlapping a TSR in the specified sample will not be associated with any TSR.
#'
#' TSS samples can be associated with TSR samples using a list or sample sheet.
#' 'sample_list' should be a named list of character vectors, with the names being the TSR
#'   sample names and the character vectors as the TSS samples(s) that should be associated with
#'   each TSR sample.
#' If 'use_sample_sheet' is set to true, a previously added sample sheet will be used to associate
#'   TSR samples with TSS samples.
#' The sample sheet can be added to the tsrexplorer object using the 'add_sample_sheet' function.
#' It should have 3 columns: 'replicate_id', 'tss_name', and 'tsr_name'.
#'
#'#' @seealso
#' \code{\link{add_sample_sheet}} to add a sample sheet to the tsrexplorer object.
#'
#' @rdname associate_with_tsr-function
#' @export

associate_with_tsr <- function(
  experiment,
  sample_list = NA,
  use_sample_sheet = FALSE
) {

  ## Check inputs.
  if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsrexplorer object")
  if (!is(use_sample_sheet, "logical")) stop("use_sample_sheet must be logical")
  if (
    (!use_sample_sheet && all(is.na(sample_list))) |
    (use_sample_sheet && all(!is.na(sample_list)))
  ) {
    stop("a sample list or sample sheet must be specified")
  }
  if (
    !use_sample_sheet && (!is(sample_list, "list") ||
    is.null(names(sample_list)))
  )
  {
    stop("sample_list must be a named list")
  }

  ## Grab sample sheet if being used.
  if (use_sample_sheet) {
    samples <- experiment@meta_data$sample_sheet
    sample_list <- as.list(samples$tss_name)
    sample_list <- set_names(sample_list, samples$tsr_name)
  }
  
  ## Associate TSSs with TSRs.
  associated_TSSs <- imap(sample_list, function(tss_names, tsr_name) {

    # Make GRanges of TSSs.
    tss_gr <- extract_counts(experiment, "tss", tss_names) %>%
      rbindlist(idcol = "sample") %>%
      as_granges

    # Make GRanges of TSRs.
    tsr_set <- extract_counts(experiment, "tsr", tsr_name) %>%
      rbindlist(idcol = "tsr_sample")
    
    tsr_gr <- tsr_set
    setnames(
      tsr_gr,
      old = c("FID", "FHASH", "width", "score", "n_unique"),
      new = c("TSR_FID", "TSR_FHASH", "tsr_width", "tsr_score", "tsr_n_unique")
    )
    tsr_gr <- as_granges(tsr_gr)

    ## Associate TSSs with overlapping TSRs
    overlapping <- tss_gr %>%
      join_overlap_left_directed(tsr_gr) %>%
      as.data.table %>%
      split(by="sample", keep.by=FALSE)

    return(overlapping)

  })

  ## Add TSSs back to the tsrexplorer object.
  associated_TSSs <- purrr::reduce(associated_TSSs, c)

  experiment <- set_count_slot(
    experiment, associated_TSSs,
    "counts", data_type, "raw"
  )
  return(experiment)
}
