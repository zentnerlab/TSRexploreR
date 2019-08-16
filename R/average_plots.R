
#' Average Plots
#'
#' Generate average plots of TSSs or TSRs
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr select filter between
#' @importFrom magrittr %>%
#' @importFrom purrr pmap
#'
#' @param experiment tsr_explorer object with annotated TSSs
#' @param sample Name of sample to make the average plot for
#' @param data_type Make average plot for TSS or TSR data
#' @param consider_score Should the TSS or TSR score be considered, or just the unique location.
#' @param upstream Bases upstream of plot center
#' @param downstream Bases downstream of plot center
#' @param threshold threshold value for TSSs
#'
#' @return ggplot2 object of average plot
#'
#' @export
#' @rdname plot_average-function

plot_average <- function(experiment, sample, data_type = c("tss", "tsr"), consider_score = FALSE, upstream = 1000, downstream = 1000, threshold = 1) {

	## Pull data out of appropriate slot.
	if (data_type == "tss") {
		sample_data <- experiment@annotated$TSSs[[sample]]
		color_type <- "#431352"
	} else if (data_type == "tsr") {
		sample_data <- experiment@annotated$TSRs[[sample]]
		color_type <- "#34698c"
	}

	## Preliminary preparation of data.
	sample_data <- sample_data %>%
		select(distanceToTSS, score) %>%
		filter(
			score >= threshold,
			between(distanceToTSS, -upstream, downstream)
		)

	## Grab info needed for plotting.
	if (consider_score) {
		plot_data <- sample_data %>%
			pmap(., function(distanceToTSS, score) rep(distanceToTSS, score)) %>%
			unlist %>%
			enframe(name = NULL, value = "distanceToTSS")
	} else {
		plot_data <- select(sample_data, distanceToTSS)
	}

	## Plot averages.
	p <- ggplot(annotated, aes(distanceToTSS)) +
		geom_density(fill = color_type, color = color_type) +
		labs(
			x = "TSS",
			y = "Density",
			title = sample
		) +
		theme_bw()

	return(p)
}
