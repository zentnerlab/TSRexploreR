
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
#' @param tsr_metrics Names of metrics in tsrexplorer TSR granges to plot
#' @param plot_type Output either a 'violin', 'jitter', 'box', or "boxjitter" plot (default violin)
#' @param samples Either 'all' or a vector of sample names to analyze
#' @param log2_transform Whether the metric should be log2 transformed prior to plotting
#' @param ncol Number of columns to plot data
#' @param use_cpm Whether to use the CPM normalized or raw counts
#' @param dominant Whether to only consider dominant TSRs
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
	dominant = FALSE,
	...
) {

	## Grab data.
	selected_data <- experiment %>%
		extract_counts("tsr", samples, use_cpm) %>%
		bind_rows(.id = "sample") %>%
		as.data.table

	selected_data <- selected_data[, c("sample", "dominant", tsr_metrics), with = FALSE]

	## Consider only dominant TSRs if required.
	if (dominant) {
		selected_data <- selected_data[(dominant)]
	}
	selected_data <- selected_data[, dominant := NULL]

	## Log2+1 transform data if requested
	if (log2_transform) {
		selected_data[,
			(tsr_metrics) := lapply(.SD, function(x) log2(x + 1)),
			.SDcols = tsr_metrics
		]
		setnames(
			selected_data, old = tsr_metrics,
			new = str_c("Log2(", tsr_metrics, " + 1)")
		)
		selected_data <- as.data.table(selected_data)
	}

	## Prepare data for plotting.
	if (log2_transform) {
		selected_data <- melt(
			selected_data,
			measure.vars = str_c("Log2(", tsr_metrics, " + 1)"),
			variable.name = "metric"
		)
	} else {
		selected_data <- melt(
			selected_data,
			measure.vars = tsr_metrics,
			variable.name = "metric"
		)
	}

	## Make density plot of metric.
	p <- ggplot(selected_data, aes(x = sample, y = value))

	if (plot_type == "violin") {
		p <- p + geom_violin(aes(fill = sample), ...)
	} else if (plot_type == "box") {
		p <- p + geom_boxplot(fill = NA, aes(color = sample), ...)
	} else if (plot_type == "jitter") {
		p <- p + geom_jitter(aes(color = sample), ...)
	} else if (plot_type == "boxjitter") {
		p <- p +
			geom_jitter(color = "lightgrey", ...) +
			geom_boxplot(fill = NA, aes(color = sample), outlier.shape = NA)
	}
		
	p <- p +
		theme_bw() +
		scale_fill_viridis_d() +
		scale_color_viridis_d() +
		facet_wrap(~ metric, ncol = ncol, scales = "free") +
		theme(axis.text.x = element_blank())

	return(p)
}
