#' DE MA Plot
#'
#' Generate MA plot for differential TSRs or Genes (RNA-seq)
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr case_when mutate
#'
#' @param experiment tsrexplorer object
#' @param de_comparisons Which differential expression comparisons to plot
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param ncol Number of columns for the facets
#' @param ... Arguments passed to geom_point
#'
#' @return ggplot2 object of differential TSRs volcano plot.
#'
#' @rdname plot_volcano-function
#'
#' @export

plot_ma <- function(
	experiment,
	data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	de_comparisons = "all", ncol = 1, ...
){

	## Get differential expression tables.
	de_samples <- extract_de(experiment, data_type, de_comparisons)

	## Prepare DE data for plotting.
	de_samples <- bind_rows(de_samples)
	de_samples <- de_samples[, .(sample, FID, log2FC, logCPM, DE)]
	de_samples[, DE := factor(DE, levels = c("up", "unchanged", "down"))]

	## MA plot of differential expression
	p <- ggplot(de_samples, aes(x = logCPM, y = log2FC, color = DE)) +
		geom_point(...) +
		theme_bw() +
		scale_color_viridis_d() +
		facet_wrap(~ sample, ncol = ncol, scales = "free")

	return(p)
}

#' Export to clusterProfiler
#'
#' Export DEGs for use in clusterProfiler term enrichment.
#'
#' @import tibble
#' @importFrom dplyr select mutate case_when filter
#' 
#' @param annotated_de Annotated differential TSRs
#' @param log2fc_cutoff Log2 fold change cutoff for significance
#' @param fdr_cutoff FDR cutoff for significance
#'
#' @rdname export_for_enrichment-function
#'
#' @export

export_for_enrichment <- function(annotated_de, log2fc_cutoff = 1, fdr_cutoff = 0.05) {
	
	## Prepare data for export.
	export_data <- annotated_de %>%
		select(geneId, log2FC, FDR) %>%
		mutate(change = case_when(
			log2FC >= log2fc_cutoff & FDR <= fdr_cutoff ~ "increase",
			log2FC <= -log2fc_cutoff & FDR <= fdr_cutoff ~ "decrease",
			TRUE ~ "unchanged"
		)) %>%
		filter(change != "unchanged")

	return(export_data)
}

#' Plot DE Numbers
#'
#' Plot number of DE features.
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The comparisons to plot
#' @param ... Additional arguments passed to geom_col
#'
#' @rdname plot_num_de-function
#' @export

plot_num_de <- function(
	experiment, data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	de_comparisons = "all", ...
) {
	## Get appropriate samples.
	de_samples <- experiment %>%
		extract_de(data_type, de_comparisons) %>%
		bind_rows

	## prepare data for plotting.
	de_data <- de_samples[, .(count = .N), by = .(sample, DE)]
	de_data[, DE := factor(DE, levels = c("up", "unchanged", "down"))]

	## Plot data.
	p <- ggplot(de_data, aes(x = sample, y = count, fill = DE)) +
		geom_col(position = "stack", ...) +
		theme_bw() +
		scale_fill_viridis_d() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))

	return(p)
}
