#' Merge Samples
#'
#' Merge replicates or selected samples.
#'
#' @param experiment TSRexploreR object
#' @param data_type Whether to merge TSSs or TSRs
#' @param sample_sheet Sample sheet
#' @param merge_group Column in sample sheet to merge by
#' @param merge_list Named list of samples to merge
#' @param merge_replicates If 'TRUE', replicate groups will be merged
#' @param threshold Filter out TSSs or TSRs below this raw count threshold before merging
#' @param sample_list If merge_replicates is set to 'FALSE',
#' specify what samples to merge in list format.
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

  ## If group is specified prepare sample list.
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
  samples <- map(merged_samples, as.data.table)
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
#' @param experiment TSRexploreR object
#' @param use_sample_sheet Whether to use a sample sheet as a key for association of TSS and TSR samples
#' @param sample_list If 'use_sample_sheet' is FALSE, provide a list with TSR and TSS sample names
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#' tsre_exp <- associate_with_tsr(tsre_exp, sample_list=list(
#'   "S288C_WT_1"="S288C_WT_1", "S288C_WT_2"="S288C_WT_2", "S288C_WT_3"="S288C_WT_3",
#'   "S288C_D_1"="S288C_D_1", "S288C_D_2"="S288C_D_2", "S288C_D_3"="S288C_D_3"
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
#' The sample sheet can be added to the TSRexploreR object using the 'add_sample_sheet' function.
#' It should have 3 columns: 'replicate_id', 'tss_name', and 'tsr_name'.
#'
#'#' @seealso
#' \code{\link{add_sample_sheet}} to add a sample sheet to the TSRexploreR object.
#'
#' @rdname associate_with_tsr-function
#' @export

associate_with_tsr <- function(
  experiment,
  sample_list=NA,
  use_sample_sheet=FALSE
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.flag(use_sample_sheet))

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
    sample_list <- magrittr::set_names(sample_list, samples$tsr_name)
  }
  
  ## Associate TSSs with TSRs.
  associated_TSSs <- imap(
    sample_list,
    ~ .tss_association(experiment, .x, .y)
  )

  ## Add TSSs back to the TSRexploreR object.
  associated_TSSs <- purrr::reduce(associated_TSSs, c)

  experiment <- set_count_slot(
    experiment, associated_TSSs,
    "counts", "tss", "raw"
  )
  return(experiment)
}

#' TSS Association
#'
#' @param experiment tsr explorer object
#' @param tss_names names of TSSs that will be ssociated with the TSR
#' @param tsr_name name of TSR that TSSs will be associated with

.tss_association <- function(
  experiment,
  tss_names,
  tsr_name
) {

  # Make GRanges of TSSs.
  tss_gr <- extract_counts(experiment, "tss", tss_names) %>%
    rbindlist(idcol="sample") %>%
    as_granges

  # Make GRanges of TSRs.
  tsr_set <- extract_counts(experiment, "tsr", tsr_name) %>%
    rbindlist(idcol="tsr_sample")

  tsr_gr <- tsr_set
  setnames(
    tsr_gr,
    old=c("FHASH", "width", "score", "n_unique"),
    new=c("TSR_FHASH", "tsr_width", "tsr_score", "tsr_n_unique")
  )
  tsr_gr <- as_granges(tsr_gr)

  ## Associate TSSs with overlapping TSRs
  overlapping <- tss_gr %>%
    join_overlap_left_directed(tsr_gr) %>%
    as.data.table %>%
    split(by="sample", keep.by=FALSE)

  return(overlapping)

}
