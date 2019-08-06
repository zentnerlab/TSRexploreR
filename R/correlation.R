#' Replicate Correlation Matrix
#'
#' Correlation matrix to explore replicate concordance.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#'
#' @param experiment tsrexplorer object with TMM normalized counts
#' @param corr_metric Correlation metric ("pearson", "spearman")
#'
#' @return ggplot2 object
#'
#' @export
#' @rdname correlation_matrix_plot-methods

setGeneric(
	"plot_corr_matrix",
	function(experiment, ...) standardGeneric("plot_corr_matrix")
)

#' @rdname correlation_matrix_plot-methods

setMethod(
	"plot_corr_matrix", signature="tsr_object",
	function(experiment, corr_metric=c("pearson", "spearman")) {
		## Prepare data for plotting.
		corr.matrix <- experiment@TMM %>%
			dplyr::select(-TSS_position) %>%
			as.matrix %>%
			cor(., method = corr_metric) %>%
			as_tibble(rownames = "sample_1", .name_repair="unique") %>%
			gather(key = "sample_2", value = corr_metric, -sample_1) %>%
			mutate(corr_metric = round(corr_metric, 3))
	
		## Plot correlation matrix.
		p <- ggplot(corr.matrix, aes(x=sample_1, y=sample_2, fill=corr_metric, label=corr_metric)) +
			geom_tile(color="white", lwd=0.5) +
			geom_label(color="white", label.size=NA, fill=NA) +
			scale_fill_viridis_c(limits=c(0.9,1), name=corr_metric) +
			theme_minimal() +
			theme(
				axis.text.x=element_text(angle=45, hjust=1),
				panel.grid=element_blank(),
				axis.title=element_blank()
			)

		return(p)
	}
)
