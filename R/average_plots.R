
#' TSS Average Plots
#'
#' Generate average plot of TSSs versus annotated TSSs
#'
#' @include tsrexplorer.R
#' @include annotate.R
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr select filter
#' @importFrom magrittr %>%
#' @importFrom purrr pmap
#'
#' @param experiment tsr_explorer object with annotated TSSs
#' @param upstream Bases upstream of plot center
#' @param downstream BZZases downstream of plot center
#' @param threshold threshold value for TSSs
#'
#' @return ggplot2 object of average plot
#'
#' @export
#' @rdname plot_tss_average-function

plot_tss_average <- function(experiment, sample, upstream = 1000, downstream = 1000, threshold = 1) {
	## Grab info needed for plotting.
	annotated <- experiment@annotated$TSSs[[sample]] %>%
		select(distanceToTSS, score) %>%
		filter(
			score >= threshold,
			between(distanceToTSS, -upstream, downstream)
		)

	## Split out individual TSSs from summed TSSs.
	annotated <- annotated %>%
		pmap(., function(distanceToTSS, score) rep(distanceToTSS, score)) %>%
		unlist %>%
		enframe(value = "distanceToTSS") %>%
		select(-name)

	## Plot TSS averages.
	p <- ggplot(annotated, aes(distanceToTSS)) +
		geom_density(fill = "#431352", color = "#431352") +
		labs(
			x = "Start Codon",
			y = "Density",
			title = sample
		) +
		theme_bw()

	return(p)
}

#' TSR Average Plots
#'
#' Generate average plot of TSRs versus annotated TSSs
#'
#' @include tsrexplorer.R
#' @include annotate.R
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr select filter
#' @importFrom magrittr %>%
#' @importFrom purrr pmap
#'
#' @param experiment tsr_explorer object with annotated TSRs
#' @param upstream Bases upstream of plot center
#' @param downstream Bases downstream of plot center
#'
#' @return ggplot2 object of average plot
#'
#' @export
#' @rdname plot_tsr_average-function

plot_tsr_average <- function(experiment, sample, upstream = 1000, downstream = 1000) {
	## Grab info needed for plotting.
	annotated <- experiment@annotated$TSRs[[sample]] %>%
		select(distanceToTSS, score) %>%
		filter(between(distanceToTSS, -upstream, downstream))

	## Split out individual TSRs from summed TSSs.
	annotated <- annotated %>%
		pmap(., function(distanceToTSS, score) rep(distanceToTSS, score)) %>%
		unlist %>%
		enframe(value = "distanceToTSS") %>%
		select(-name)

	## Plot TSS averages.
	p <- ggplot(annotated, aes(distanceToTSS)) +
		geom_density(color = "#34698c", fill = "#34698c") +
		labs(
			x = "Start Codon",
			y = "Density",
			title = sample
		) +
	theme_bw()

	return(p)
}
