#'
#' heatmaps and/or scatter plots to explore replicate concordance of TSSs or TSRs.
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate_at left_join
#' @importFrom purrr discard
#' @importFrom tidyr gather
#' @importFrom GGally ggpairs
#'
#' @param experiment tsrexplorer object with TMM normalized counts
#' @param data_type Whether to make scatter plots from TSS or TSR data
#' @param correlation_plot Whether to make a correlation 'heatmap', 'scatter', or 'combined'
#' @param correlation_metric Use either spearman or pearson correlation
#' @param samples Either "all" or vector of samples to plot
#' @param log2 Should the TMM values be log2+1 transformed prior to plotting?
#'
#' @return ggplot2 object
#'
#' @export
#' @rdname plot_correlation-function

plot_correlation <- function(experiment, data_type = c("tss", "tsr", "rnaseq_v_tss"), correlation_plot = "combined", correlation_metric = "pearson", samples = "all", log2_transform = TRUE) {
	
	## Get data from proper slot.
	if (data_type == "tss") {
		normalized_counts <- experiment@normalized_counts$TSSs
		type_color <- "#431352"
	} else if (data_type == "tsr") {
		normalized_counts <- experiment@normalized_counts$TSRs
		type_color <- "#34698c"
	} else if (data_type == "rnaseq_v_tss") {
		normalized_counts <- left_join(
			experiment@normalized_counts$RNAseq_features,
			experiment@normalized_counts$TSS_features,
			by = "position"
		)
		type_color <- "#29AF7FFF"
	}

	## Log2+1 transform data if indicated.
	pre_transformed <- normalized_counts
	if (log2_transform) normalized_counts <- mutate_at(normalized_counts, vars(-position), ~log2(. + 1))

	## Select all samples if "all" specified.
	if (samples == "all") {
		samples <- normalized_counts %>%
			colnames %>%
			discard(. == "position")
	}

	## Make functions for scatter plot.

	# Create custom scatter plot format.
	custom_scatter <- function(data, mapping) {
		ggplot(data = data, mapping = mapping) +
			geom_point(size = 0.25, color = type_color) +
			geom_abline(intercept = 0, slope = 1, lty = 2)
	}

	# Create custom correlation heatmap format.
	custom_heatmap <- function(data, mapping) {
		sample_1 <- pre_transformed[ ,as.character(mapping$x)[2]]
		sample_2 <- pre_transformed[ ,as.character(mapping$y)[2]]

		correlation <- cor(sample_1, sample_2, method = correlation_metric) %>%
			round(3) %>%
			as_tibble(name_repair = "unique", rownames = "sample_1") %>%
			gather(key = "sample_2", value = "correlation", -sample_1)

		ggplot(correlation, aes(x = sample_1, y = sample_2)) +
			geom_tile(color = "white", aes(fill = correlation)) +
			geom_label(aes(label = correlation), label.size=NA, fill=NA, color = "white") +
			scale_fill_viridis_c(limits = c(0, 1), direction = -1)
	}

	## Plot the correlation plot. 	

	if (correlation_plot == "scatter") {
		p <- ggpairs(
			normalized_counts,
			columns = samples,
			upper = list(continuous = custom_scatter),
			lower = NULL,
			diag = NULL
		)
	} else if (correlation_plot == "heatmap") {
		p <- ggpairs(
			normalized_counts,
			columns = samples,
			upper = list(continuous = custom_heatmap),
			lower = NULL,
			diag = NULL
		)
	} else if (correlation_plot == "combined") {
		p <- ggpairs(
			normalized_counts,
			columns = samples,
			upper = list(continuous = custom_heatmap),
			lower = list(continuous = custom_scatter)
		)
	}

	return(p)
}
