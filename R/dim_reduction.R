
#' Dimensionality Reduction.
#'
#' @description
#' Dimensionality reduction plot using PCA or UMAP.
#'
#' @importFrom uwot umap
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param method 'umap' or 'pca'
#' @param n_neighbors Nearest neighbors parameter for umap
#' @param min_dist Minimum distance between points parameter for umap
#' @param ... Additional arguments passed to geom_point
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
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#' plot_reduction(tsre_exp, data_type = "tsr")
#'
#' @return ggplot2 plot of dimension reduction
#'
#' @rdname plot_reduction-function
#' @export

plot_reduction <- function(
	experiment,
	data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	method = "umap",
	n_neighbors = 2,
	min_dist = 0.1,
	...
) {

	## Input checks.
	if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsrexplorer object")

        if (!is(data_type, "character") || length(data_type) > 1) {
                stop("data_type must be a character")
        }
        data_type <- str_to_lower(data_type)
        if (!data_type %in% c("tss", "tsr", "tss_features", "tsr_features")) {
		stop("data_type must be 'tss', 'tsr', 'tss_features', or 'tsr_features'")
	}

	if (!is(method, "character") || length(method) > 1) stop("method must be 'umap' or 'pca'")
	method <- str_to_lower(method)
	if (!method %in% c("umap", "pca")) stop("method must be 'umap' or 'pca'")

	if (
		!is(n_neighbors, "numeric") || length(n_neighbors) > 1 ||
		n_neighbors %% 1 != 0 || n_neighbors < 2
	) {
		stop("n_neighbors must be a positive integer greater than or equal to 2")
	}

	if (
		!is(min_dist, "numeric") || length(min_dist) > 1 || min_dist < 0
	){
		stop("min_dist must be a number greater than or equal to 0")
	}

	## Grab TMM counts.
	tmm_counts <- extract_matrix(experiment, data_type, "all")

	## Filter out counts with too little expression.
	keep_index <- tmm_counts %>%
		assay("counts") %>%
		filterByExpr

	tmm_counts <- tmm_counts[keep_index, ]

	## Prepare data for umap.
	tmm_counts <- tmm_counts %>%
		assay("tmm") %>%
		t %>%
		as.data.frame %>%
		rownames_to_column("sample") %>%
		as.data.table

	## UMAP dimensionality reduction.
	umap_results <- tmm_counts %>%
		umap(n_neighbors = n_neighbors, min_dist = min_dist)
	
	umap_results <- as.data.table(umap_results)
	setnames(umap_results, old = c("V1", "V2"), new = c("UMAP_1", "UMAP_2"))
	umap_results[, sample := tmm_counts[["sample"]]]

	## Plot UMAP.
	p <- ggplot(umap_results, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
		geom_point(...) +
		theme_bw() +
		scale_color_viridis_d()

	return(p)
}
