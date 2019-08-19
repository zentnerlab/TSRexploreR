
#' Average Plots
#'
#' Generate average plots of TSSs or TSRs
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr select filter between bind_rows
#' @importFrom magrittr %>%
#' @importFrom purrr pmap map
#'
#' @param experiment tsr_explorer object with annotated TSSs
#' @param samples Either 'all' to plot all samples, or a vector of sample names
#' @param data_type Make average plot for TSS or TSR data
#' @param consider_score Should the TSS or TSR score be considered, or just the unique location.
#' @param upstream Bases upstream of plot center
#' @param downstream Bases downstream of plot center
#' @param threshold threshold value for TSSs
#' @param ncol Number of columns to use when plotting data
#'
#' @return ggplot2 object of average plot
#'
#' @export
#' @rdname plot_average-function

plot_average <- function(experiment, data_type = c("tss", "tsr"), samples = "all", consider_score = FALSE, upstream = 1000, downstream = 1000, threshold = 1, ncol = 1) {

	## Pull data out of appropriate slot.
	if (data_type == "tss") {
		if (samples == "all") samples <- names(experiment@annotated$TSSs)
		sample_data <- experiment@annotated$TSSs[samples] %>%
			bind_rows(.id = "samples")
		color_type <- "#431352"
	} else if (data_type == "tsr") {
		if (samples == "all") samples <- names(experiment@annotated$TSRs)
		sample_data <- experiment@annotated$TSRs[samples] %>%
			bind_rows(.id = "samples") %>%
			rename(score = nTAGs)
		color_type <- "#34698c"
	}

	## Preliminary preparation of data.
	sample_data <- sample_data %>%
		select(distanceToTSS, score, samples) %>%
		filter(
			score >= threshold,
			between(distanceToTSS, -upstream, downstream)
		)

	## Grab info needed for plotting.
	if (consider_score) {
		plot_data <- sample_data %>%
			split(.$samples) %>%
			map(., ~pmap(., function(distanceToTSS, score, samples) rep(distanceToTSS, score)) %>%
					unlist %>%
					enframe(name = NULL, value = "distanceToTSS")
			) %>%
			bind_rows(.id = "samples")
	} else {
		plot_data <- select(sample_data, samples, distanceToTSS)
	}

	## Plot averages.
	p <- ggplot(plot_data, aes(distanceToTSS)) +
		geom_density(fill = color_type, color = color_type) +
		labs(
			x = "TSS",
			y = "Density"
		) +
		theme_bw() +
		facet_wrap(~ samples)

	return(p)
}
