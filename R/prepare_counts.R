#' Format Counts
#'
#' @description
#' Format TSS or TSR counts for further analysis. 
#'
#' @inheritParams common_params
#' @param data_type Whether to format TSS or TSR counts 
#'
#' @details
#' When TSSs or TSRs are first loaded into the TSRexploreR object
#'   they are stored as GRanges objects.
#' This function converts these into data.table format,
#'   and adds a few important columns for downstream analysis.
#'
#' @return TSRexploreR object with properly formatted features
#'    in data.table format.
#'
#' @examples
#' data(TSSs)
#'
#' TSSs[1] %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss")
#'
#' @export

format_counts <- function(
  experiment,
  data_type=c("tss", "tsr"),
  samples="all"
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr"))
  assert_that(is.character(samples))

  ## Get selected samples and generate raw count matrices.
  if (data_type == "tss") {

    ## Get selected samples.
    select_samples <- tss_experiment(experiment)
    if (!all(samples == "all")) select_samples <- select_samples[samples]

  } else if (data_type == "tsr") {

    ## Get selected samples.
    select_samples <- tsr_experiment(experiment)
    if (!all(samples == "all")) select_samples <- select_samples[samples]

  }

  ## Turn counts into data.table
  if (data_type %in% c("tss", "tsr")) {
    raw_counts <- map(select_samples, function(x) {
      x <- as.data.table(x)
      x[, FHASH := str_c(seqnames, start, end, strand, sep=":")]
      return(x)
    })
  }

  ## Place counts in proper TSRexploreR object slot.
  if (data_type == "tss") {
    experiment@counts$TSSs$raw <- c(experiment@counts$TSSs$raw, raw_counts)
  } else if (data_type == "tsr") {
    experiment@counts$TSRs$raw <- c(experiment@counts$TSRs$raw, raw_counts)
  }

  return(experiment)
}
