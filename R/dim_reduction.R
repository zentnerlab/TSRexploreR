
#' Dimensionality Reduction.
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
#' @rdname plot_reduction-function
#' @export

plot_reduction <- function(
	experiment, data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	method = "umap", n_neighbors = 2, min_dist = 0.1, ...
) {
	
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
