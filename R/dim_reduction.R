#' Dimensionality Reduction.
#'
#' @description
#' Dimensionality reduction plot using PCA
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#' @param remove_var Remove features in this bottom fraction.
#' @param center Center the data.
#' @param scale Scale the data.
#' @param ... Additional arguments passed to PCAtools::biplot.
#'
#' @return ggplot2 object of PCA plot.
#'
#' @details
#' This function will generatete a PCA plot of the first two PCs.
#' This helps to visualize the relative similarity of samples
#'   based on the most variables features.
#'
#' 'remove_var' removes features in the bottom fraction of variance.
#' 'center' and 'scale' will center and scale the data respectively.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#'
#' tsre_exp <- TSSs %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss")
#'
#' # TSS PCA plot.
#' \donttest{plot_reduction(tsre, data_type="tss")}
#'
#' # TSR PCA plot.
#' tsre <- tss_clustering(tsre, threshold=3)
#' \donttest{plot_reduction(tsre, data_type="tsr")}
#'
#' @export

plot_reduction <- function(
  experiment,
  samples="all",
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  use_normalized=TRUE,
  remove_var=NULL,
  center=TRUE,
  scale=TRUE,
  ...
) {

  ## Check if PCAtools is installed.
  if (!requireNamespace("PCAtools", quietly = TRUE)) {
    stop("Package \"PCAtools\" needed for this function to work. Please install it.",
      call. = FALSE)
  }

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))
  assert_that(is.flag(use_normalized))
  assert_that(is.null(remove_var) || (remove_var > 0 & remove_var < 1))

  ## Grab counts.
  count_mat <- experiment %>%
    extract_counts(data_type, samples) %>%
    .count_matrix(data_type, use_normalized)

  ## Grab sample sheet
  metadata <- copy(experiment@meta_data$sample_sheet)
  metadata[, c("file_1", "file_2") := NULL]
  metadata <- metadata[sample_name %in% colnames(count_mat)]
  metadata <- metadata[match(metadata[, sample_name], colnames(count_mat))]
  metadata <- column_to_rownames(metadata, "sample_name")

  ## Create biplot.
  p <- count_mat %>%
    PCAtools::pca(
      center=center, scale=scale, metadata=metadata,
      removeVar=remove_var
    ) %>%
    PCAtools::biplot(...)

  return(p) 
}
