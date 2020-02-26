
#' Mark Dominant
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter group_by ungroup between bind_rows
#' @importFrom purrr map
#'
#' @param experiment tsrexplorer object with annotated TSSs/TSRs
#' @param data_type Either 'tss' or 'tsr'
#' @param threshold Read threshold for TSS/TSRs
#' @param mark_per By default marks dominant TSR per gene, and dominant TSS per TSR.
#' TSSs can also be set per as dominant TSS per 'gene'.
#'
#' @return Tibble with dominant TSSs for each gene
#'
#' @rdname mark_dominant-function
#'
#' @export

mark_dominant <- function(experiment, data_type = c("tss", "tsr"), threshold = 1, mark_per = "default") {
	
	## Select samples.
	select_samples <- extract_counts(experiment, data_type, "all")

	## Mark dominant TSS/TSR per gene if requested.
	if (data_type == "tsr" | (data_type == "tss" & mark_per == "gene")) {
		dominant <- map(select_samples, function(x) {
			x[,
				dominant := (
					score == max(score) &
					!simple_annotations %in% c("Downstream", "Intergenic") &
					score >= threshold
				),
				by = eval(ifelse(
					experiment@settings$annotation[, feature_type] == "transcript",
					"transcriptId", "geneId"
				))
			]

			return(x)
		})
	
	## Mark the dominant TSS per TSR if requested.
	} else if (data_type == "tss" & mark_per == "default") {
		dominant <- map(select_samples, function(x) {
			x[,
				dominant := (
					!is.na(score) &
					score == max(score) &
					score >= threshold
				),
				by = TSR_FID
			]

			return(x)
		})
	}

	## Return dominant TSS/TSR.
	if (data_type == "tss") {
		experiment@counts$TSSs$raw <- dominant
	} else if (data_type == "tsr") {
		experiment@counts$TSRs$raw <- dominant
	}

	return(experiment)
}

#' Max UTR Length
#'
#' Get TSS with furthest distance
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom dplyr filter between group_by select summarize bind_rows
#' @importFrom purrr map
#'
#' @param experiment tsrexplorer object with annotated TSSs
#' @param samples Either 'all' or names of sample to analyze
#' @param threshold Number of reads required for each TSS
#' @param max_upstream Max upstream distance of TSS to consider
#' @param max_downstream Max downstream distance of TSS to consider
#' @param feature_type Feature type used when finding distance to TSS ("geneId", "transcriptId")
#' @param quantiles Number of quantiles to break data into.
#'
#' @return tibble with max UTR length for features
#'
#' @rdname max_utr-function
#'
#' @export

max_utr <- function(
	experiment,
	samples = "all",
	threshold = 1,
	max_upstream = 1000,
	max_downstream = 100,
	feature_type = c("geneId", "transcriptId"),
	quantiles = NA
) {
	## Grab selected samples.
	max_utr <- experiment %>%
		extract_counts("tss", samples) %>%
		bind_rows(.id = "sample") %>%
		as.data.table

	setnames(max_utr, old = feature_type, new = "feature")

	## Filter data.
	max_utr <- max_utr[
		score >= threshold &
		dplyr::between(distanceToTSS, -max_upstream, max_downstream),
		.(sample, distanceToTSS, feature, score)
	]

	## Get TSS with minimum distance to start codon.
	max_utr <- max_utr[,
		.SD[which.min(distanceToTSS)],
		by = .(sample, feature)
	]

	## Add quantiles if requested.
	if (!is.na(quantiles)) {
		max_utr[, ntile := ntile(score, quantiles), by = sample]
	}

	## Return DataFrame.
	max_utr <- DataFrame(max_utr)
	metadata(max_utr)$quantiles <- quantiles
	metadata(max_utr)$max_upstream <- max_upstream
	metadata(max_utr)$max_downstream <- max_downstream
	metadata(max_utr)$threshold <- threshold

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
#' @param ncol Number of columns to plot the data to
#' @param ... Arguments passed to geom_density
#'
#' @return ggplot2 object of max UTR length average
#'
#' @rdname plot_max_utr-function
#'
#' @export

plot_max_utr <- function(
	max_utr, upstream = 1000, downstream = 100,
	ncol = 1, consider_score = FALSE, ...
) {

	## Grab some info from DataFrame.
	quantiles <- metadata(max_utr)$quantiles

	## Prepare data for plotting.
	utr_plot <- as.data.table(max_utr)
	
	if (!is.na(quantiles)) {
		utr_plot[, ntile := fct_rev(factor(ntile))]
	}

	## Format data if score should be considered.
	if (consider_score) {
		utr_plot <- utr_plot[rep(seq_length(.N), score)]
	}

	## Plot max UTR length detected.
	p <- ggplot(utr_plot, aes(x = distanceToTSS)) +
		geom_density(fill = "#431352", color = "#431352")+#, ...) +
		xlim(-upstream, downstream) +
		theme_bw() +
		labs(
			x = "Max UTR Length",
			y = "Density"
		) +
		geom_vline(xintercept = 0, lty = 2)

	if (!is.na(quantiles)) {
		p <- p + facet_grid(ntile ~ sample)
	} else {
		p <- p + facet_wrap(~ sample, ncol = ncol)
	}

	return(p)
}
