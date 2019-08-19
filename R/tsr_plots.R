
#' Plot TSR Stats
#'
#' Plot selected TSR stats
#'
#' @import tibble
#' @import ggplot2
#' @importFrom purrr pluck
#' @importFrom dplyr mutate_at select
#'
#' @param experiment tsrexplorer object with TSR granges
#' @param sample Name of sample to return plot for
#' @param metric Name of metric in tsrexplorer TSR granges to plot
#' @param log2_transform Whether the metric should be log2+1 transformed prior to plotting
#' @param xlims Vector of values specifying the upper and lower limits of the x-axis
#'
#' @return ggplot2 object with tsr metrix plotted
#'
#' @export
#' @rdname plot_tsr_metric-function

plot_tsr_metric <- function(experiment, sample, metric, log2_transform = FALSE, xlims = c()) {

	## Grab data.
	selected_data <- experiment@experiment$TSRs %>%
		pluck(sample) %>%
		as_tibble(.name_repair = "unique") %>%
		select(metric)

	x_label <- metric

	## Log2+1 transform data if requested
	if (log2_transform) {
		selected_data <- mutate_at(selected_data, metric, ~log2(. + 1))
		x_label <- paste0("Log2(", metric, " + 1)")
	}

	## Make density plot of metric.
	p <- ggplot(selected_data, aes_string(metric)) +
		geom_density(fill = "#34698c", color = "#34698c") +
		theme_bw() +
		labs(
			x = x_label,
			y = "Density"
		)

	## Add xlim if set.
	if (!is.null(xlims)) {
		p <- p + xlim(xlims)
	}

	return(p)
}
