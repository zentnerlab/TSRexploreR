
#' Threshold Exploration
#'
#' Explore various naive thresholds
#'
#' @import tibble
#' @importFrom purrr map_df
#'
#' @param experiment tsrexplorer object with TSSs
#' @param samples Either 'all' or the names of the samples to analyze
#' @param max_threshold Threshold from 1 to max_threshold will be explored
#'
#' @return tibble containing information for each threshold and sample
#'
#' @rdname explore_thresholds-function
#'
#' @export

explore_thresholds <- function(
	experiment,
	max_threshold = 50,
	samples = "all"
) {
	## Grab settings information.
	feature_type <- experiment@settings$annotation[["feature_type"]]
	feature_type <- ifelse(feature_type == "transcript", "transcriptId", "geneId")

	## Pull out appropriate samples.
	select_samples <- extract_counts(experiment, "tss", samples)
	select_samples <- rbindlist(select_samples, idcol = "sample")

	## Grab information needed for thresholding plot.
	summarized_data <- map_df(seq_len(max_threshold), function(x) {
		filtered <- select_samples[score >= x]

		feature_stats <- filtered[,
			.(count = uniqueN(get(feature_type))),
			by = .(sample, promoter_proximity = ifelse(
				simple_annotations == "Promoter",
				"n_promoter_proximal", "n_promoter_distal"
			))
		]

		feature_stats <- dcast(
			feature_stats, sample ~ promoter_proximity,
			value.var = "count"
		)

		feature_stats[,
			c("n_total", "frac_promoter_proximal", "threshold") := list(
				n_promoter_proximal + n_promoter_distal,
				n_promoter_proximal / (n_promoter_proximal + n_promoter_distal),
				x
			)
		]
	})

	return(summarized_data)
}

#' Plot Threshold Exploration
#'
#' Make a plot to explore threshold values
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate
#'
#' @param threshold_data Tibble of threshold exploration data from explore_thresholds
#' @param ncol Number of columns to plot data
#' @param point_size The size of the points on the plot
#' @param ... Arguments passed to geom_point
#'
#' @return ggplot2 object containing the threshold exploration plot
#'
#' @rdname plot_threshold_exploration-function
#'
#' @export

plot_threshold_exploration <- function(threshold_data, ncol = 1, point_size = 1, ...) {

	## Plot data.
	p <- ggplot(threshold_data, aes(x = threshold, y = frac_promoter_proximal)) +
		geom_line(color = "lightgrey") +
		geom_point(aes(color = n_total), size = point_size, ...) +
		scale_color_viridis_c() +
		theme_bw() +
		xlab("Count Threshold") +
		ylab("Fraction of Promoter Proximal TSSs") +
		facet_wrap(. ~ sample, ncol = ncol)

	return(p)
}
