
#' Dominant TSS Average Distance
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter group_by ungroup between bind_rows
#'
#' @param experiment tsrexplorer object with annotated TSSs
#' @param samples Either 'all' or name of samples to analyze
#' @param threshold Read threshold for TSS
#' @param max_upstream Max upstream distance of TSS to consider
#' @param max_downstream Max downstream distance of TSS to consider
#' @param feature_type Feature that was used to find TSS distance ("geneId", "transcriptId")
#'
#' @return Tibble with dominant TSSs for each gene
#'
#' @rdname dominant_tss-function
#'
#' @export

dominant_tss <- function(
	experiment,
	samples = "all",
	threshold = 1, 
	max_upstream = 2000,
	max_downstream = 500, 
	feature_type = c("transcriptId", "geneId")
) {
	## Select samples.
	if (samples == "all") samples <- names(experiment@annotated$TSSs$raw)
	select_samples <- experiment@annotated$TSSs$raw[samples]	

	## Get dominant TSSs.
	dominant <- select_samples %>%
		map(
			~ rename(.x, "feature" = feature_type) %>%
				select(feature, distanceToTSS, score) %>%
				filter(
					score >= threshold,
					between(distanceToTSS, -max_upstream, max_downstream)
				) %>%
				group_by(feature) %>%
				filter(score == max(score)) %>%
				ungroup
		) %>%
		bind_rows(.id = "sample")

	return(dominant)
}

#' Plot Dominant TSS
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr count filter select
#' @importFrom purrr pmap
#'
#' @param dominant_tss Tibble of dominant TSSs from dominant_tss
#' @param upstream Bases upstream to display on plot
#' @param downstream Bases downstream to display on plot
#' @param ncol Number of columns to plot the data to
#'
#' @return ggplot2 object dominant TSS plot
#'
#' @rdname plot_dominant_tss-function
#'
#' @export

plot_dominant_tss <- function(dominant_tss, upstream = 2000, downstream = 500, ncol = 1) {
	## Format data for plotting
	#dominant_tss <- dominant_tss %>%
		#count(distanceToTSS) %>%
		#pmap(function(distanceToTSS, n) rep(distanceToTSS, n)) %>%
		#unlist %>%
		#enframe(value="distanceToTSS") %>%
		#select(-name)

	## Plot data
	p <- ggplot(dominant_tss, aes(distanceToTSS)) +
		geom_density(fill="#431352", color="#431352") +
		xlim(-upstream, downstream) +
		theme_bw() +
		labs(
			x="Start Codon",
			y="Density"
		) +
		geom_vline(xintercept=0, lty=2) +
		facet_wrap(~ sample, ncol = ncol)

	return(p)
}

#' Max UTR Length
#'
#' Get TSS with furthest distance
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom dplyr filter between group_by_at select summarize
#'
#' @param experiment tsrexplorer object with annotated TSSs
#' @param sample Name of sample to analyze
#' @param threshold Number of reads required for each TSS
#' @param max_upstream Max upstream distance of TSS to consider
#' @param max_downstream Max downstream distance of TSS to consider
#' @param feature_type Feature type used when finding distance to TSS ("geneId", "transcriptId")
#'
#' @return tibble with max UTR length for features
#'
#' @rdname max_utr-function
#'
#' @export

max_utr <- function(
	experiment, sample, threshold = 1,
	max_upstream = 1000, max_downstream = 100,
	feature_type = c("geneId", "transcriptId")
) {
	## Get TSS with minimum distance to start codon.
	max_utr <- experiment@annotated$TSSs[[sample]] %>%
		select(feature_type, distanceToTSS, score) %>%
		filter(
			score >= threshold,
			between(distanceToTSS, -max_upstream, max_downstream)
		) %>%
		group_by_at(feature_type) %>%
		summarize(tss_max_distance = min(distanceToTSS))

	return(max_utr)
}

#' Plot Max UTR Length
#'
#' Plot TSS with furthest distance
#'
#' @import tibble
#' @import ggplot2
#'
#' @param max_utr tibble of max UTRs output by max_utr
#' @param upstream Bases upstream to extend average to
#' @param downstream Bases downstream to extend average to
#'
#' @return ggplot2 object of max UTR length average
#'
#' @rdname plot_max_utr-function
#'
#' @export

plot_max_utr <- function(max_utr, upstream = 1000, downstream = 100) {
	p <- ggplot(max_utr, aes(x = tss_max_distance)) +
		geom_density(fill = "#431352", color = "#431352") +
		xlim(-upstream, downstream) +
		theme_bw() +
		labs(
			x = "Max UTR Length",
			y = "Density"
		) +
		geom_vline(xintercept = 0, lty = 2)

	return(p)
}
