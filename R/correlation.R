#' Replicate TSS Correlation Matrix
#'
#' Correlation matrix to explore replicate concordance.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr select mutate
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TMM normalized counts
#' @param corr_metric Correlation metric ("pearson", "spearman")
#'
#' @return ggplot2 object
#'
#' @export
#' @rdname plot_tss_corr-function

plot_tss_corr <- function(experiment, corr_metric=c("pearson", "spearman")) {
	## Prepare data for plotting.
	corr.matrix <- experiment@normalized_counts$TSSs %>%
		select(-TSS_position) %>%
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

#' Replicate TSS Scatter Plot
#'
#' Scatter plots to explore replicate concordance.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import ggplot2
#'
#' @param experiment tsrexplorer object with TMM normalized counts
#' @param sample_1 First sample name to plot
#' @param sample_2 Second sample name to plot
#'
#' @return ggplot2 object
#'
#' @export
#' @rdname plot_tss_scatter-function

plot_tss_scatter <- function(experiment, sample_1, sample_2) {
	p <- ggplot(experiment@normalized_counts$TSSs, aes_string(x=sample_1, y=sample_2)) +
		geom_point(size=0.25, color="#431352") +
		theme_bw() +
		scale_fill_viridis_d() +
		geom_abline(intercept=0, slope=1, lty=2)

	return(p)
}

#' Replicate TSR Correlation Matrix
#'
#' Correlation matrix to explore replicate concordance.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr select mutate
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TMM normalized TSR counts
#' @param corr_metric Correlation metric ("pearson", "spearman")
#'
#' @return ggplot2 object
#'
#' @export
#' @rdname plot_tsr_corr-function

plot_tsr_corr <- function(experiment, corr_metric=c("pearson", "spearman")) {
	## Prepare data for plotting.
	corr.matrix <- experiment@normalized_counts$TSRs %>%
		select(-consensus_TSR) %>%
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

#' Replicate TSR Scatter Plot
#'
#' Scatter plots to explore replicate concordance.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import ggplot2
#'
#' @param experiment tsrexplorer object with TMM normalized TSR counts
#' @param sample_1 First sample name to plot
#' @param sample_2 Second sample name to plot
#'
#' @return ggplot2 object
#'
#' @export
#' @rdname plot_tsr_scatter-function

plot_tsr_scatter <- function(experiment, sample_1, sample_2) {
	## Log2 transform value
	log2_counts <- experiment@normalized_counts$TSRs %>%
		mutate_if(is.numeric, ~log2(. + 1))

	## Plot TSR correlation scatter plot
	p <- ggplot(log2_counts, aes_string(x=sample_1, y=sample_2)) +
		geom_point(size=0.25, color = "#34698c", fill = "#34698c") +
		theme_bw() +
		scale_fill_viridis_d() +
		geom_abline(intercept=0, slope=1, lty=2)

	return(p)
}
