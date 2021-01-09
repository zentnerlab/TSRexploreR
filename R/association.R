#' Merge Samples
#'
#' Merge replicates or selected samples.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#' @param merge_group Column in sample sheet to merge by.
#' @param merge_list Named list of samples to merge.
#' @param merge_replicates If 'TRUE', replicate groups will be merged.
#' @param sample_list If merge_replicates is set to 'FALSE',
#'   specify what samples to merge in list format.
#'
#' @rdname merge_samples-function
#' @export

merge_samples <- function(
  experiment,
  data_type=c("tss", "tsr"),
  threshold=NULL,
  sample_sheet=NULL,
  merge_group=NULL,
  merge_list=NULL
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr")
  )
  assert_that(
    is.null(threshold) ||
    (is.numeric(threshold) && threshold >= 0)
  )
  assert_that(
    is.null(sample_sheet) ||
    (is.character(sample_sheet) | is.data.frame(sample_sheet))
  )
  assert_that(is.null(merge_group) || is.string(merge_group))
  assert_that(
    is.null(merge_list) ||
    (is.list(merge_list) && has_attr(merge_list, "names"))
  )
  if (is.null(merge_group) & is.null(merge_list)) {
    stop("Either 'merge_group' or 'merge_list' must be specified")
  }
  if (!is.null(merge_group) & !is.null(merge_list)) {
    stop("Either 'merge_group' or 'merge_list' must be specified")
  }

  ## If group is specified, prepare sample list.
  merge_type <- case_when(
    !is.null(merge_group) ~ "sample_sheet",
    !is.null(merge_list) ~ "sample_list"
  )

  if (merge_type == "sample_sheet") {
    keep <- c("sample_name", merge_group)
    merge_list <- experiment@meta_data$sample_sheet[, ..keep] %>%
      split(by=merge_group, keep.by=FALSE) %>%
      map(`[[`, "sample_name")
  }

  ## Merge feature sets.
  merged_samples <- map(merge_list, function(x) {

    # Bind ranges.
    samples <- experiment %>%
      extract_counts(data_type, x) %>%
      rbindlist %>%
      as_granges

    # Merge overlapping ranges.
    if ("normalized_score" %in% colnames(elementMetadata(samples))) {
      samples <- reduce_ranges_directed(
        samples, score=sum(score),
        normalized_score=sum(normalized_score)
      )
    } else {
      samples <- reduce_ranges_directed(samples, score=sum(score))
    }

    # Sort ranges.
    samples <- sort(samples)

    return(samples)
    
  })

  ## Add merged GRanges to experiment slot.
  if (data_type == "tss") {
    experiment@experiment$TSSs <- c(experiment@experiment$TSSs, merged_samples)
  } else if (data_type == "tsr") {
    experiment@experiment$TSRs <- c(experiment@experiment$TSRs, merged_samples)
  }

  ## Add merged GRanges to count slot.
  merged_samples <- map(merged_samples, as.data.table)
  if (data_type == "tss") {
    experiment@counts$TSSs$raw <- c(experiment@counts$TSSs$raw, merged_samples)
  } else if (data_type == "tsr") {
    experiment@counts$TSRs$raw <- c(experiment@counts$TSRs$raw, merged_samples)
  }

  ## Add sample info to sample sheet.
  if (is.null(experiment@meta_data$sample_sheet)) {
    setDT(sample_sheet)
    experiment@meta_data$sample_sheet <- sample_sheet
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
#' @inheritParams common_params
#' @param sample_list List with TSRs as names and TSSs as vector.
#'   If NULL will associate TSSs with TSRs from the sample of the same name.
#'
#' @details
#' TSRexploreR provides many options for the import and merging of TSSs and TSRs.
#' Because of this, TSS samples must be associated with TSR samples after TSR import,
#'   TSR merging, or TSS clustering using this function.
#' This adds an extra workflow step, but provides more analytical flexibility.
#' Each TSS with genomic coordinates overlapping those of a TSR in the specified TSR sample
#'   will be linked to that TSR.
#' TSSs not overlapping a TSR in the specified sample will not be associated with any TSR.
#'
#' TSS samples can be associated with TSR samples using a list or sample sheet.
#' 'sample_list' should be a named list of character vectors, with the names being the TSR
#'   sample names and the character vectors as the TSS samples(s) that should be associated with
#'   each TSR sample. If no sample list is provided, the function will automatically associate
#'   TSSs with the TSRs from the sample of the same name.
#'   
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#' exp <- format_counts(exp, data_type="tss")
#' exp <- tss_clustering(exp)
#' exp <- associate_with_tsr(exp, sample_list=list(
#'   "S288C_WT_1"="S288C_WT_1", "S288C_WT_2"="S288C_WT_2", "S288C_WT_3"="S288C_WT_3",
#'   "S288C_D_1"="S288C_D_1", "S288C_D_2"="S288C_D_2", "S288C_D_3"="S288C_D_3"
#' ))
#'
#'#' @seealso
#' \code{\link{add_sample_sheet}} to add a sample sheet to the TSRexploreR object.
#'
#' @rdname associate_with_tsr-function
#' @export

associate_with_tsr <- function(
  experiment,
  sample_list=NULL
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is.null(sample_list) ||
    (is.list(sample_list) && has_attr(sample_list, "names"))
  )

  ## Get TSSs.
  if (!is.null(sample_list)) {
    tss <- map(sample_list, function(samples) {
      map(samples, ~extract_counts(experiment, "tss", .x, FALSE))
    })
  } else {
    tss <- experiment@counts$TSSs$raw %>%
      imap(~purrr::set_names(list(.x), .y))
  }

  ## Get TSRs.
  if (!is.null(sample_list)) {
    tsr <- map(
      sample_list,
      ~extract_counts(experiment, "tsr", .x, FALSE)
    )
  } else {
    tsr <- copy(experiment@counts$TSRs$raw)
  }

  ## Associate TSSs with TSRs.
  tss <- imap(tsr, function(tsr, tsr_name) {
    tsr[, tsr_sample := tsr_name]
    setkey(tsr, seqnames, strand, start, end)
    tss <- tss[[tsr_name]]
    tss <- map(tss, function(x) {
      setkey(x, seqnames, strand, start, end)
      overlap <- foverlaps(x, tsr)
      overlap[, c("start", "end") := NULL]
      setnames(
        overlap,
        old=c(
          "width", "n_unique", "FHASH", "i.start", "i.end",
          "i.width", "i.FHASH", "score", "i.score"
        ),
        new=c(
          "tsr_width", "tsr_n_unique", "TSR_FHASH", "start",
          "end", "width", "FHASH", "tsr_score", "score"
        )
      )
      if (any(colnames(x) == "normalized_score")) {
        setnames(
          overlap,
          old=c("normalized_score", "i.normalized_score"),
          new=c("tsr_normalized_score", "normalized_score")
        )
      }
      overlap <- overlap %>%
        as_granges %>%
        sort %>%
        as.data.table
      return(overlap)
    })
  })
  tss <- flatten(tss)

  ## Add TSSs back to the TSRexploreR object.
  experiment <- set_count_slot(
    experiment, tss, "counts", "tss", "raw"
  )
  return(experiment)
}
