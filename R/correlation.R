
#' Generate Correlation Matrix
#'
#' Calculate a correlation matrix for sample concordance of TSSs or TSRs
#'
#' @import tibble
#' @importFrom purrr discard
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom dplyr select_at vars
#' @importFrom tidyr gather
#'
#' @param experiment tsrexplorer object with TMM normalized counts
#' @param data_type Whether to make scatter plots from TSS, TSR, or RNA-seq & five-prime data
#' @param samples Either 'all' or the names of the samples to plot
#' @param correlation_metric Use either spearman or pearson correlation
#'
#' @return correlation matrix
#'
#' @rdname find_correlation-function
#'
#' @export

find_correlation <- function(
	experiment,
	data_type = c("tss", "tsr", "features"),
	samples = "all",
	correlation_metric = "pearson"
) {
	## Select appropriate data.
	if (data_type == "tss") {
		normalized_counts <- experiment@counts$TSSs$normalized %>% assay("tmm")
		type_color <- "#431352"
	} else if (data_type == "tsr") {
		normalized_counts <- experiment@counts$TSRs$normalized %>% assay("tmm")
		type_color <- "#34698c"
	} else if (data_type == "features") {
		normalized_counts <- experiment@counts$features$normalized %>% assay("tmm")
		type_color <- "#29AF7FFF"
	}

	## Select all samples if "all" specified.
	if (length(samples) == 1 & samples == "all") samples <- colnames(normalized_counts)

	## make correlation matrix.
	correlation <- normalized_counts %>%
		.[, samples] %>%
		cor(method = correlation_metric) %>%
		as_tibble(.name_repair = "unique", rownames = "sample_1") %>%
		gather(-sample_1, key = "sample_2", value = "cor")

	## Place correlation values into proper slot.
	if (data_type == "tss") {
		experiment@correlation$TSSs <- correlation
	} else if (data_type == "tsr") {
		experiment@correlation$TSRs <- correlation
	} else {
		experiment@correlation$features <- correlation
	}

	return(experiment)		
}

#' Plot Sample Correlation
#'
#' heatmaps and/or scatter plots to explore replicate concordance of TSSs or TSRs.
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate_at left_join vars
#' @importFrom purrr discard
#' @importFrom tidyr gather
#' @importFrom GGally ggpairs
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2 
#' @importFrom viridis viridis
#' @importFrom grid gpar
#'
#' @param experiment tsrexplorer object with TMM normalized counts
#' @param data_type Whether to make scatter plots from TSS, TSR, or RNA-seq & five-prime data
#' @param samples Either 'all' or the names of the samples to plot
#' @param correlation_plot Whether to make a correlation 'heatmap', 'scatter', 'combined' or 'hierarchical'
#' @param correlation_metric Use either spearman or pearson correlation
#' @param log2_transform Should the TMM values be log2+1 transformed prior to plotting?
#' @param font_size The font size for the heatmap tiles
#' @param pt_size Point size for the scatter plots
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap or GGally::ggpairs
#'
#' @return ggplot2 object
#'
#' @rdname plot_correlation-function
#'
#' @export

plot_correlation <- function(
	experiment,
	data_type = c("tss", "tsr", "features"),
	samples = "all",
	correlation_plot = "combined",
	correlation_metric = "pearson",
	log2_transform = TRUE,
	font_size = 4,
	pt_size = 0.5,
	...
) {
	
	## Get data from proper slot.
	if (data_type == "tss") {
		normalized_counts <- experiment@counts$TSSs$tmm_matrix
		type_color <- "#431352"
	} else if (data_type == "tsr") {
		normalized_counts <- experiment@counts$TSRs$tmm_matrix
		type_color <- "#34698c"
	} else if (data_type == "features") {
		normalized_counts <- experiment@counts$features$tmm_matrix
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
			geom_point(size = pt_size, color = type_color, stroke = 0) +
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
			geom_label(aes(label = correlation), label.size=NA, fill=NA, color = "black", size = font_size) +
			scale_fill_viridis_c(limits = c(0, 1))
	}

	## Plot the correlation plot. 	

	if (correlation_plot == "scatter") {
		p <- ggpairs(
			normalized_counts,
			columns = samples,
			upper = list(continuous = custom_scatter),
			lower = NULL,
			diag = NULL,
			...
		)
	} else if (correlation_plot == "heatmap") {
		p <- ggpairs(
			normalized_counts,
			columns = samples,
			upper = list(continuous = custom_heatmap),
			lower = NULL,
			diag = NULL,
			...
		)
	} else if (correlation_plot == "combined") {
		p <- ggpairs(
			normalized_counts,
			columns = samples,
			upper = list(continuous = custom_heatmap),
			lower = list(continuous = custom_scatter),
			...
		)
	} else if (correlation_plot == "hierarchical") {
		corr_matrix <- pre_transformed %>%
			column_to_rownames("position") %>%
			cor(method = correlation_metric)

		p <- Heatmap(
			corr_matrix,
			name = correlation_metric,
			row_names_gp = gpar(fontsize = font_size),
			column_names_gp = gpar(fontsize = font_size),
			heatmap_legend_param = list(
				title_gp = gpar(fontsize = font_size),
				labels_gp = gpar(fontsize = font_size),
				grid_height = unit(2, "mm"),
				grid_width = unit(3, "mm")
			),
			...
		)
	}

	return(p)
}
