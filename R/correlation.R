#' Replicate TSS Correlation Matrix
#'
#' Correlation matrix to explore replicate concordance.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr select mutate rename
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TMM normalized counts
#' @param data_type Whether to create a correlation matrix for TSSs, TSRs, or RNA-seq
#' @param corr_metric Correlation metric ("pearson", "spearman")
#'
#' @return ggplot2 object
#'
#' @export
#' @rdname plot_correlation-function

plot_correlation <- function(experiment, data_type = c("tss", "tsr", "rnaseq"), corr_metric=c("pearson", "spearman")) {

	## Select data from tsrexplorer object.
	if (data_type == "tss") {
		corr_matrix <- experiment@normalized_counts$TSSs
	} else if (data_type == "tsr") {
		corr_matrix <- experiment@normalized_counts$TSRs
	} else if (data_type == "rnaseq") {
		corr_matrix <- experiment@normalized_counts$RNAseq
	}

	## Prepare data for plotting.
	corr_matrix <- corr_matrix %>%
		rename("feature" = 1) %>%
		column_to_rownames("feature") %>%
		as.matrix %>%
		cor(method = corr_metric) %>%
		as_tibble(rownames = "sample_1", .name_repair="unique") %>%
		gather(key = "sample_2", value = corr_metric, -sample_1) %>%
		mutate(corr_metric = round(corr_metric, 3))
	
	## Plot correlation matrix.
	p <- ggplot(corr_matrix, aes(x=sample_1, y=sample_2, fill=corr_metric, label=corr_metric)) +
		geom_tile(color="white", lwd=0.5) +
		geom_label(color="white", label.size=NA, fill=NA) +
		scale_fill_viridis_c(limits=c(0,1), name=corr_metric, direction = -1) +
		theme_minimal() +
		theme(
			axis.text.x=element_text(angle=45, hjust=1),
			panel.grid=element_blank(),
			axis.title=element_blank()
		)

	return(p)
}

#' Replicate Scatter Plots
#'
#' Scatter plots to explore replicate concordance of TSSs or TSRs.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import ggplot2
#'
#' @param experiment tsrexplorer object with TMM normalized counts
#' @param sample_1 First sample name to plot
#' @param sample_2 Second sample name to plot
#' @param data_type Whether to make scatter plots from TSS or TSR data
#'
#' @return ggplot2 object
#'
#' @export
#' @rdname plot_scatter-function

plot_scatter <- function(experiment, sample_1, sample_2, data_type = c("tss", "tsr")) {
	
	## Get data from proper slot.
	if (data_type == "tss") {
		normalized_counts <- experiment@normalized_counts$TSSs
		type_color = "#431352"
	} else if (data_type == "tsr") {
		normalized_counts <- experiment@normalized_counts$TSRs
		type_color = "#34698c"
	}

	p <- ggplot(normalized_counts, aes_string(x = sample_1, y = sample_2)) +
		geom_point(size = 0.25, color = type_color) +
		theme_bw() +
		scale_fill_viridis_d() +
		geom_abline(intercept = 0, slope = 1, lty = 2)

	return(p)
}
