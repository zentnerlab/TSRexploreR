#' CPM Normalize Counts
#'
#' @description
#' CPM normalize the TSS, TSR, and/or feature counts.
#'
#' @param experiment tsrexplorer object
#' @param data_type 'tss', 'tsr', 'tss_features', or 'tsr_features'
#'
#' @details
#' Counts Per Million (CPM) is an approach for read number normalization, 
#'   allowing comparison between samples for various downstream analyses and plots.
#' For more quantitative comparisons, TMM normalization is recommended (this should be fleshed out more qq)
#'
#' @return tsrexplorer object with CPM values.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- cpm_normalize(tsre_exp, data_type = "tss")
#'
#' @rdname cpm_normalize-function
#' @export

cpm_normalize <- function(
  experiment,
  data_type = c("tss", "tsr", "tss_features", "tsr_features")
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))

  ## Get selected samples.
  select_samples <- switch(
    data_type,
    "tss"=experiment@counts$TSSs$raw,
    "tsr"=experiment@counts$TSRs$raw,
    "tss_features"=experiment@counts$TSS_features$raw,
    "tsr_features"=experiment@counts$TSR_features$raw
  )

  ## CPM normalize counts.
  cpm_counts <- select_samples %>%
    map(function(x) {
      x[, cpm := cpm(score)]
      return(x)
    })

  ## Add CPM-normalized counts back to tsrexplorer object.
  experiment <- set_count_slot(
    experiment, cpm_counts,
    "counts", data_type, "raw"
  )

  return(experiment)
}

#' TMM Normalize TSSs or TSRs
#'
#' @description
#' Using edgeR to TMM normalize TSSs or TSRs.
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#'
#' @param experiment tsrexplorer object
#' @param data_type Whether TSSs, TSRs, or RNA-seq & 5' feature counts should be normalized
#' @param threshold Consider only features with at least this number of raw counts
#' @param n_samples Filter out positions without features meeting the the selected threshold
#'   in this number of samples
#' @param samples Vector with names of samples to include in the normalization
#'
#' @details
#' The TMM normalization method, employed by edgeR, is designed to reduce the influence of
#'   library size on inter-sample comparisons.
#'
#' For TMM normalization, it is recommended to remove features with few or no reads as
#'   these may bias the final results.
#' To facilitate this filtering, two arguments are provided, 'threshold' and 'n_samples'.
#' Features must have greater than or equal to 'threshold' number of raw counts in at least
#'   'n_samples' number of samples to be retained.
#'
#' @return tsrexplorer object with tmm normalized count matrices
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- count_matrix(tsre_exp, data_type = "tss")
#' tsre_exp <- tmm_normalize(exp, data_type = "tss")
#'
#' @seealso \code{\link{count_matrix}} to prepare the matrices.
#'   \code{\link{plot_correlation}} for various correlation plots.
#'
#' @rdname tmm_normalize-function
#' @export

tmm_normalize <- function(
  experiment,
  data_type = c("tss", "tsr", "tss_features", "tsr_features"),
  samples = "all",
  threshold = 1,
  n_samples = 1
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features")
  )
  assert_that(is.character(samples))
  assert_that(is.count(threshold))
  assert_that(is.count(n_samples))

  ## Get selected samples.
  select_samples <- extract_matrix(experiment, data_type, samples)

  ## Filter counts.
  sample_matrix <- as.data.table(assay(select_samples, "counts"))
  sample_matrix[, match := rowSums((.SD >= threshold)) >= n_samples]
  keep_ids <- which(sample_matrix[, match])
  
  select_samples <- select_samples[keep_ids, ]
  filtered_counts <- assay(select_samples, "counts")

  ## TMM normalize filtered counts.
  tmm_counts <- filtered_counts %>%
    DGEList %>%
    calcNormFactors %>%
    cpm

  ## Create filtered and TMM normalized RangedSummarizedExperiment.
  assay(select_samples, "tmm") <- tmm_counts

  ## Add TMM-normalized counts to tsrexplorer object.
  experiment <- set_count_slot(
    experiment, select_samples,
    "counts", data_type, "matrix"
  )

  return(experiment)
}
