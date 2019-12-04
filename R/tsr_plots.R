
#' Plot TSR Stats
#'
#' Plot selected TSR stats
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate_at select select_at bind_rows vars mutate
#' @importFrom tidyr gather
#'
#' @param experiment tsrexplorer object with TSR granges
#' @param sample Name of sample to return plot for
#' @param tsr_metrics Names of metrics in tsrexplorer TSR granges to plot
#' @param plot_type Output either a 'violin', 'jitter', 'box', or "boxjitter" plot (default violin)
#' @param samples Either 'all' or a vector of sample names to analyze
#' @param log2_transform Whether the metric should be log2 transformed prior to plotting
#' @param ncol Number of columns to plot data
#' @param use_cpm Whether to use the CPM normalized or raw counts
#' @param ... Arguments passed to ggplot2 plotting functions
#'
#' @return ggplot2 object with tsr metrix plotted
#'
#' @rdname plot_tsr_metric-function
#'
#' @export

plot_tsr_metric <- function(
	experiment,
	tsr_metrics,
	plot_type = "violin",
	samples = "all",
	log2_transform = FALSE,
	ncol = 1,
	use_cpm = FALSE,
	...
) {

	## Grab data.
	if (use_cpm) {
		if (samples == "all") samples <- names(experiment@annotated$TSRs$cpm)
		select_samples <- experiment@annotated$TSRs$cpm[samples]
	} else {
		if (samples == "all") samples <- names(experiment@annotated$TSRs$raw)
		select_samples <- experiment@annotated$TSRs$raw[samples]
	}

	selected_data <- select_samples %>%
		bind_rows(.id = "samples") %>%
		select_at(vars("samples", tsr_metrics)) %>%
		gather(key = "tsr_stat", value = "stat_value", -samples)

	## Log2+1 transform data if requested
	if (log2_transform) {
		selected_data <- selected_data %>%
			mutate(
				stat_value = log2(stat_value),
				tsr_stat = paste0("Log2(", tsr_stat, ")")
			)
	}

	## Make density plot of metric.
	p <- ggplot(selected_data, aes(x = samples, y = stat_value))

	if (plot_type == "violin") {
		p <- p + geom_violin(aes(fill = samples), ...)
	} else if (plot_type == "box") {
		p <- p + geom_boxplot(fill = NA, aes(color = samples), ...)
	} else if (plot_type == "jitter") {
		p <- p + geom_jitter(aes(color = samples), ...)
	} else if (plot_type == "boxjitter") {
		p <- p +
			geom_jitter(color = "lightgrey", ...) +
			geom_boxplot(fill = NA, aes(color = samples), outlier.shape = NA)
	}
		
	p <- p +
		theme_bw() +
		scale_fill_viridis_d() +
		scale_color_viridis_d() +
		facet_wrap(~ tsr_stat, ncol = ncol, scales = "free") +
		theme(axis.text.x = element_blank())

	return(p)
}
