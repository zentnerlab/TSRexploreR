
#' Threshold Exploration
#'
#' Explore various naive thresholds.
#'
#' @importFrom purrr map_df
#'
#' @param experiment tsrexplorer object with TSSs
#' @param samples Either 'all' or the names of the samples to analyze
#' @param max_threshold Threshold from 1 to max_threshold will be explored
#'
#' @return data.frame containing information for each threshold and sample
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data = annotation,
#'   data_type = "tss", feature_type = "transcript"
#' )
#' thresh_results <- explore_thresholds(tsre_exp)
#'
#' @seealso \code{\link{plot_threshold_exploration}} to plot the results.
#'
#' @rdname explore_thresholds-function
#' @export

explore_thresholds <- function(
	experiment,
	max_threshold = 50,
	samples = "all"
) {
	## Check inputs.
	if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")

	if (!is(max_threshold, "numeric")) stop("max_threshold must be a positive integer")
	if (max_threshold %% 1 != 0) stop("max_threshold must be a positive integer")
	if (max_threshold < 5) stop("max_threshold must be greater than or equal to 5")

	if (!is(samples, "character")) stop("samples must be a character vector")

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
#' Make a plot to explore threshold values.
#'
#' @import ggplot2
#'
#' @param threshold_data Tibble of threshold exploration data from explore_thresholds
#' @param ncol Number of columns to plot data
#' @param point_size The size of the points on the plot
#' @param ... Arguments passed to geom_point
#'
#' @return ggplot2 object containing the threshold exploration plot
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data = annotation,
#'   data_type = "tss", feature_type = "transcript"
#' )
#' thresh_results <- explore_thresholds(tsre_exp)
#' plot_threshold_exploration(thresh_results)
#'
#' @seealso \code{\link{explore_thresholds}} for initial calculations.
#'
#' @rdname plot_threshold_exploration-function
#' @export

plot_threshold_exploration <- function(threshold_data, ncol = 1, point_size = 1, ...) {

	## Check inputs.
	if (!is(threshold_data, "data.frame")) stop("threshold_data must be a data.frame")
	
	if (!is(ncol, "numeric")) stop("ncol must be a positive integer")
	if (ncol %% 1 != 0) stop("ncol must be a positive integer")
	if (ncol < 1) stop("ncol must be a positive integer")

	if (!is(point_size, "numeric")) stop("point_size must be a positive number")
	if (!point_size > 0) stop("point_size must be a positive number")

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
