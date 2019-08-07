
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
#' @importFrom purrr map pmap
#'
#' @param experiment tsr_explorer object with annotated TSSs
#' @param upstream Bases upstream of plot center
#' @param downstream BZZases downstream of plot center
#' @param threshold threshold value for TSSs
#'
#' @return ggplot2 object of average plot
#'
#' @export
#' @rdname plot_tss_average

plot_tss_average <- function(experiment, sample, upstream = 1000, downstream = 1000, threshold = 1) {
	## Grab info needed for plotting.
	annotated <- map(
		experiment@annotated$TSSs,
		~select(., distanceToTSS, score) %>%
			filter(score >= threshold)
	)

	## Split out individual TSSs from summed TSSs.
	annotated <- map(
		annotated,
		~pmap(., function(distanceToTSS, score) rep(distanceToTSS, score)) %>%
			unlist %>%
			enframe(value = "distanceToTSS") %>%
			select(-name)
	)

	## Plot TSS averages.
	p <- ggplot(annotated[[sample]], aes(distanceToTSS))
		geom_density(fill = "#431352", color = "#431352") +
		xlim(-upstream, downstream) +
		labs(
			x = "Start Codon",
			y = "Density",
			title = sample
		) +
		theme_bw()

	return(p)
}
