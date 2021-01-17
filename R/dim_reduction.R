#' Dimensionality Reduction.
#'
#' @description
#' Dimensionality reduction plot using PCA
#'
#' @importFrom PCAtools pca biplot
#'
#' @inheritParams common_params
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param remove_var Remove features in this bottom fraction
#' @param center center the data
#' @param scale scale the data
#' @param ... Additional arguments passed to PCAtools::biplot
#'
#' @details
#' This function will generate a dimensionality reduction plot.
#' These help to visualize the relative similarity of samples
#'   based on the most variables features.
#'
#' 'method' lets you choose between the more traditional PCA plot or newer
#'   UMAP plot.
#'
#' If a UMAP plot is chosen, there are two parameters available
#'   for tweaking the results.
#' 'n_neighbors' will find the specified number of nearest neighbor points,
#'   and 'min_dist' specifies the minimum distance allowed between two points 
#'   in the UMAP plot. 
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#' plot_reduction(tsre_exp, data_type="tsr")
#'
#' @return ggplot2 plot of dimension reduction
#'
#' @rdname plot_reduction-function
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

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))
  assert_that(is.flag(use_normalized))
  assert_that(is.null(remove_var) | (remove_var > 0 & remove_var < 1))

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
    pca(
      center=center, scale=scale, metadata=metadata,
      removeVar=remove_var
    ) %>%
    biplot(...)

  return(p) 
}
