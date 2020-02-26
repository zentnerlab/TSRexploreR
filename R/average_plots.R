
#' Average Plots
#'
#' Generate average plots of TSSs or TSRs
#'
#' @import tibble
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr bind_rows bind_cols select filter between group_by mutate ntile ungroup
#' @importFrom magrittr %>% extract
#' @importFrom purrr map pmap
#' @importFrom forcats fct_rev
#' @importFrom tidyr separate
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowRanges
#'
#' @param experiment tsr_explorer object with annotated TSSs
#' @param samples Either 'all' to plot all samples, or a vector of sample names
#' @param data_type Make average plot for TSS or TSR data
#' @param consider_score Should the TSS or TSR score be considered, or just the unique location.
#' @param upstream Bases upstream of plot center
#' @param downstream Bases downstream of plot center
#' @param threshold threshold value for TSSs
#' @param ncol Number of columns to use when plotting data when quantiles not set
#' @param quantiles Number of quantiles to split data into
#' @param use_cpm Whether to use the CPM normalized or raw counts
#' @param dominant Consider only dominant TSS or TSR
#' @param ... Arguments passed to geom_density
#'
#' @return ggplot2 object of average plot
#'
#' @rdname plot_average-function
#'
#' @export

plot_average <- function(
	experiment,
	data_type = c("tss", "tsr"),
	samples = "all",
	consider_score = FALSE,
	upstream = 1000,
	downstream = 1000,
	threshold = 1,
	ncol = 1,
	quantiles = NA,
	use_cpm = FALSE,
	dominant = FALSE,
	...
) {

	## Assign color type.
	if (data_type == "tss") {
		color_type <- "#431352"
	} else if (data_type == "tsr") {
		color_type <- "#34698c"
	}

	## Pull data out of appropriate slot.
	sample_data <- experiment %>%
		extract_counts(data_type, samples, use_cpm) %>%
		bind_rows(.id = "samples") %>%
		as.data.table

	keep_cols <- c("distanceToTSS", "score", "samples")
	if (dominant) keep_cols <- c(keep_cols, "dominant")

	sample_data <- sample_data[
		score >= threshold &
		dplyr::between(distanceToTSS, -upstream, downstream),
		..keep_cols, 
	]

	## Only keep dominant if requested.
	if (dominant) {
		sample_data <- sample_data[(dominant)]
		sample_data[, dominant := NULL]
	}

	## Add quantile info if requested.
	if (!is.na(quantiles)) {
		sample_data[,
			ntile := ntile(score, quantiles),
			by = .(samples, distanceToTSS)
		]
	}

	## Update data if score is also considered instead of just unique position.
	if (consider_score) sample_data <- sample_data[rep(seq_len(.N), score)]

	## Plot averages.
	p <- ggplot(sample_data, aes(distanceToTSS)) +
		geom_density(fill = color_type, color = color_type, ...) +
		labs(
			x = "Position Relative to Annotated TSS",
			y = "Density"
		) +
		theme_bw()

	if (!is.na(quantiles)) {
		p <- p + facet_grid(fct_rev(factor(ntile)) ~ samples)
	} else {
		p <- p + facet_wrap(~ samples, ncol = ncol)
	}
	return(p)
}
