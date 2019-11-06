
#' Average Plots
#'
#' Generate average plots of TSSs or TSRs
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr select filter between bind_rows ntile group_by ungroup
#' @importFrom magrittr %>%
#' @importFrom purrr pmap map
#' @importFrom forcats fct_rev
#'
#' @param experiment tsr_explorer object with annotated TSSs
#' @param samples Either 'all' to plot all samples, or a vector of sample names
#' @param data_type Make average plot for TSS or TSR data
#' @param consider_score Should the TSS or TSR score be considered, or just the unique location.
#' @param upstream Bases upstream of plot center
#' @param downstream Bases downstream of plot center
#' @param threshold threshold value for TSSs
#' @param ncol Number of columns to use when plotting data when quantiles not set
#' @param quantiles Number of quantiles to split data into
#'
#' @return ggplot2 object of average plot
#'
#' @rdname plot_average-function
#'
#' @export

plot_average <- function(
	experiment,
	data_type = c("tss", "tsr"),
	samples = "all",
	consider_score = FALSE,
	upstream = 1000,
	downstream = 1000,
	threshold = 1,
	ncol = 1,
	quantiles = 1
) {

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

	## Add quantile info if requested.
	sample_data <- sample_data %>%
		group_by(samples, distanceToTSS) %>%
		mutate(ntile = ntile(score, quantiles)) %>%
		ungroup

	## Grab info needed for plotting.
	if (consider_score) {
		sample_data <- sample_data %>%
			mutate(samples_ntile = paste0(samples, "_:_", ntile)) %>%
			select(-samples, -ntile) %>%
			split(.$samples_ntile) %>%
			map(., ~pmap(., function(distanceToTSS, score, samples_ntile) rep(distanceToTSS, score)) %>%
					unlist %>%
					enframe(name = NULL, value = "distanceToTSS")
			) %>%
			bind_rows(.id = "samples") %>%
			separate(samples, into = c("samples", "ntile"), sep = "_:_")
	}

	## Plot averages.
	p <- ggplot(sample_data, aes(distanceToTSS)) +
		geom_density(fill = color_type, color = color_type) +
		labs(
			x = "Position Relative to Annotated TSS",
			y = "Density"
		) +
		theme_bw()

	if (quantiles > 1) {
		p <- p + facet_grid(fct_rev(factor(ntile)) ~ samples)
	} else {
		p <- p + facet_wrap(~ samples, ncol = ncol)
	}
	return(p)
}
