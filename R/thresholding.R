
#' Threshold Exploration
#'
#' Explore various naive thresholds
#'
#' @import tibble
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom ChIPseeker annotatePeak
#' @importFrom dplyr group_by ungroup mutate bind_rows filter pull n_distinct
#' @importFrom purrr map
#'
#' @param experiment tsrexplorer object with TSSs
#' @param annotation_file Either path and file name of annotation, or a TxDb object
#' @param feature_type Annotate on gene or transcript level
#' @param samples Either 'all' or the names of the samples to analyze
#' @param max_threshold Threshold from 1 to max_threshold will be explored
#' @param upstream Bases upstream of TSS
#' @param downstream Bases downstream of TSS
#'
#' @return tibble containing information for each threshold and sample
#'
#' @rdname explore_thresholds-function
#'
#' @export

explore_thresholds <- function(
	experiment,
	annotation_file,
	feature_type = c("gene", "transcript"),
	max_threshold = 50,
	upstream = 1000,
	downstream = 100,
	samples = "all"
) {
	## Pull out appropriate samples.
	if (all(samples == "all")) {
		samples <- names(experiment@experiment$TSSs)
	}

	select_samples <- experiment@experiment$TSSs[samples]

	## Prepare annotation file.
	if (class(annotation_file) == "character") {
		annotation <- makeTxDbFromGFF(annotation_file)
	} else if (class(annotation_file) == "TxDb") {
		annotation <- annotation_file
	}

	## Change 'nTAGs' to 'score' if present.
	select_samples <- select_samples %>%
		map(function(x) {
			metacol_names <- x %>% mcols %>% names

			if (any(metacol_names == "nTAGs")) {
				x$score <- x$nTAGs
			}

			return(x)
		})

	## Annotate TSSs.
	raw_annotated <- select_samples %>%
		map(
			~ annotatePeak(.x,
				tssRegion = c(-upstream, downstream),
				TxDb = annotation,
				sameStrand = TRUE,
				level = feature_type,		
			) %>% as_tibble(.name_repair = "unique")
		)

	## Grab information needed for thresholding plot.
	summarize_data <- function(dataset, max) {
		set_stats <- map(1:max, function(x) {
			filtered_dataset <- filter(dataset, score >= x)

			total_genes <- filtered_dataset %>%
				pull(geneId) %>%
				n_distinct

			genes_with_promoter_proximal_tss <- filtered_dataset %>%
				filter(annotation == "Promoter") %>%
				pull(geneId) %>%
				n_distinct

			total_TSSs <- nrow(filtered_dataset)

			promoter_proximal_TSSs <- filtered_dataset %>%
				filter(annotation == "Promoter") %>%
				nrow

			frac_promoter_proximal_TSSs <- promoter_proximal_TSSs / total_TSSs

			result <- tibble(
				"threshold" = x,
				"total_genes" = total_genes,
				"genes_with_promoter_proximal_tss" = genes_with_promoter_proximal_tss,
				"total_TSSs" = total_TSSs,
				"promoter_proximal_TSSs" = promoter_proximal_TSSs,
				"frac_promoter_proximal_TSSs" = frac_promoter_proximal_TSSs
			)

			return(result)
		}) %>% bind_rows
		
		return(set_stats)
	}

	threshold_data <- raw_annotated %>%
		map(~ summarize_data(.x, max_threshold)) %>%
		bind_rows(.id = "sample")
		
	return(threshold_data)
}

#' Plot Threshold Exploration
#'
#' Make a plot to explore threshold values
#'
#' @import tibble
#' @import ggplot2
#'
#' @param threshold_data Tibble of threshold exploration data from explore_thresholds
#' @param ncol Number of columns to plot data
#' @param point_size The size of the points on the plot
#' @param sample_order Optional vector specifying sample order in plot
#' @param ... Arguments passed to geom_point
#'
#' @return ggplot2 object containing the threshold exploration plot
#'
#' @rdname plot_threshold_exploration-function
#'
#' @export

plot_threshold_exploration <- function(threshold_data, ncol = 1, point_size = 1, sample_order = NULL, ...) {

	## Change sample order if specified.
	if (!(is.null(sample_order))) {
		threshold_data <- mutate(threshold_data, sample = factor(sample, sample_order))
	}

	## Plot data.
	p <- ggplot(threshold_data, aes(x = threshold, y = frac_promoter_proximal_TSSs)) +
		geom_line(color = "lightgrey") +
		geom_point(aes(color = genes_with_promoter_proximal_tss), size = point_size, ...) +
		scale_color_viridis_c() +
		theme_bw() +
		facet_wrap(. ~ sample, ncol = ncol)

	return(p)
}
