#' Plot Sample Correlation
#'
#' Heatmaps and/or scatter plots to explore replicate concordance of TSSs or TSRs.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2 
#' @importFrom viridis viridis
#' @importFrom grid gpar grid.text
#'
#' @inheritParams common_params
#' @param data_type Whether to make scatter plots from TSS or TSRs.
#' @param correlation_metric Whether to use Spearman or Pearson correlation.
#' @param font_size The font size for the heatmap tiles.
#' @param cluster_samples Logical for whether hierarchical clustering should be performed
#'   on rows and columns.
#' @param heatmap_colors Vector of colors for heatmap.
#' @param show_values Logical for whether to show correlation values on the heatmap.
#' @param return_matrix Return the correlation matrix without plotting correlation heatmap.
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap.
#'
#' @details
#' Correlation plots are a good way to assess the similarity of samples.
#' This can be useful to determine replicate concordance and for the initial assessment of
#'   differences between samples from different conditions.
#' This function generates various plots using a previously TMM- or MOR-normalized count matrix.
#
#' Pearson correlation is recommended for samples from the same technology due to 
#' the expectation of a roughly linear relationship between the magnitudes of values 
#' for each feature. Spearman correlation is recommended for samples from different technologies,
#' such as STRIPE-seq vs. CAGE, due to the expectation of a roughly linear relationship 
#' between the ranks, rather than the specific values, of each feature.
#'
#' @return ggplot2 object of correlation heatmap,
#'   or correlation matrix if 'return_matrix' is TRUE.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss") %>%
#'   normalize_counts(data_type="tss", method="CPM")
#'
#' # TSS correlation.
#' \donttest{plot_correlation(tsre, data_type="tss")}
#'
#' # TSR correlation.
#' tsre <- tss_clustering(tsre, threshold=3)
#' \donttest{plot_correlation(tsre, data_type="tsr")}
#'
#' @seealso \code{\link{normalize_counts}} for TSS and TSR normalization.
#'
#' @export

plot_correlation <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  samples="all",
  correlation_metric="pearson",
  threshold=NULL,
  use_normalized=TRUE,
  font_size=12,
  cluster_samples=FALSE,
  heatmap_colors=NULL,
  show_values=TRUE,
  return_matrix=FALSE,
  ...
) {

  ## Check inputs.
  if (!is(experiment, "tsr_explorer")) stop("experiment must be a TSRexploreR object")
  data_type <- match.arg(data_type, c("tss", "tsr", "tss_features", "tsr_features"))
  assert_that(is.character(samples))
  correlation_metric <- match.arg(
   str_to_lower(correlation_metric),
    c("pearson", "spearman")
  )
  assert_that(is.numeric(font_size) && font_size > 0)
  assert_that(is.flag(cluster_samples))
  assert_that(is.null(heatmap_colors) | is.character(heatmap_colors))
  assert_that(is.flag(show_values))
  assert_that(is.flag(use_normalized))
  assert_that(
    is.null(threshold) ||
    (is.numeric(threshold) && threshold > 0)
  )
  assert_that(is.flag(return_matrix))

  ## Get data from proper slot.
  normalized_counts <- experiment %>%
    extract_counts(data_type, samples) %>%
    .count_matrix("tss", use_normalized)
  
  sample_names <- colnames(normalized_counts)

  ## Define default color palette.
  color_palette <- switch(
    data_type,
    "tss"="#431352",
    "tsr"="#34698c",
    "tss_features"="#29AF7FFF",
    "tsr_features"="#29AF7FFF"
  )

#  ## Log2 + 1 transform data if indicated.
#  pre_transformed <- copy(normalized_counts)
#  if (log2_transform) {
#    normalized_counts <- log2(normalized_counts + 1)
#  }

  ## Correlation Matrix.
  cor_mat <- cor(normalized_counts, method=correlation_metric)

  ## Return correlation matrix if requested.
  if (return_matrix) return(cor_mat)

  ## ComplexHeatmap Correlation Plot.
  heatmap_args <- list(
    cor_mat,
    row_names_gp=gpar(fontsize=font_size),
    column_names_gp=gpar(fontsize=font_size)
  )
  if (!cluster_samples) {
    heatmap_args <- c(heatmap_args, list(cluster_rows=FALSE, cluster_columns=FALSE))
  }
  if (!is.null(heatmap_colors)) {
    heatmap_args <- c(heatmap_args, list(col=heatmap_colors))
  }
  if (show_values) {
    heatmap_args <- c(heatmap_args, list(
      cell_fun=function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp=gpar(fontsize=font_size))
      }
    ))
  }

  p <- do.call(Heatmap, heatmap_args)

  return(p)

}
