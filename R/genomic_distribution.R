
#' Genomic Distribution
#'
#' Get genomic distribution of TSSs and TSRs
#'
#' @import tibble
#' @import data.table
#' @importFrom GenomicRanges GRanges score "score<-"
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame "metadata<-"
#' @importFrom dplyr bind_rows select filter group_by mutate ntile case_when count ungroup arrange rename bind_cols
#' @importFrom magrittr %>% extract
#' @importFrom purrr map
#' @importFrom tidyr spread
#'
#' @param experiment tsrexplorer object with annotated TSRs
#' @param samples Either "all" or vector of sample names
#' @param data_type Whether to get distribution of TSSs or TSRs
#' @param threshold Filter out TSSs or TSRs under a certain read count number
#' @param quantiles number of quantiles to break the data down into
#' @param dominant Should dominant TSS per 'gene' or 'tsr' only be considered
#'
#' @return tibble with TSS or TSR genomic distribution stats
#'
#' @rdname genomic_distribution-function
#'
#' @export 

genomic_distribution <- function(
	experiment, data_type = c("tss", "tsr"), samples = "all",
	threshold = 1, quantiles = NA, dominant = FALSE
) {

	## Extract samples.
	selected_samples <- experiment %>%
		extract_counts(data_type, samples) %>%
		bind_rows(.id = "samples")

	## Retain dominant
	keep_cols <- c("samples", "score", "simple_annotations")
	if (dominant) keep_cols <- c(keep_cols, "dominant")
 
	selected_samples <- selected_samples[
		score >= threshold,
		..keep_cols
	]
	selected_samples[,
		simple_annotations := factor(
			simple_annotations,
			levels = c("Promoter", "Exon", "Intron", "Downstream", "Intergenic")
		)
	]

	## Break data into quantiles if quantiles set is greater than 1.
	if (!is.na(quantiles)) {
		selected_samples[, ntile := ntile(score, quantiles), by = samples]
	}
	
	## Prepare data to be plotted later.
	if (!is.na(quantiles)) {
		genomic_distribution <- selected_samples[,
			.(count = .N),
			by = .(samples, simple_annotations, ntile)
		][,
			.(simple_annotations, count, fraction = count / sum(count)),
			by = .(samples, ntile)
		]
	} else {
		genomic_distribution <- selected_samples[,
			.(count = .N),
			by = .(samples, simple_annotations)
		][,
			.(simple_annotations, count, fraction = count / sum(count)),
			by = samples
		]
	}

	## Prepare DataFrame to return.
	dist_exp <- DataFrame(genomic_distribution)

	## Add quantile information to summarized experiment.
	metadata(dist_exp)$quantiles <- quantiles
	metadata(dist_exp)$dominant <- dominant

	return(dist_exp)
}

#' Plot Genomic Distribution
#'
#' Plot genomic distribution of TSSs or TSRs
#'
#' @import tibble
#' @import ggplot2
#' @importFrom forcats fct_rev
#' @importFrom dplyr mutate
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom S4Vectors metadata
#'
#' @param genomic_distribution tibble of TSS or TSR genomic distributions from tsr_genomic_distribution
#' @param sample_order Optional vector specifying order of samples to plot by sample name
#'
#' @return ggplot2 object with TSS or TSR genomic distribution plot
#'
#' @rdname plot_genomic_distribution-function
#'
#' @export 

plot_genomic_distribution <- function(genomic_distribution, sample_order = NA) {
	
	## Pull out information from DataFrame.
	genomic_dist <- as_tibble(genomic_distribution, .name_repair = "unique")

	## Order the samples if required.
	if (!is.na(sample_order)) {
		genomic_dist <- genomic_dist %>%
			mutate(samples = fct_rev(factor(samples, levels = sample_order)))
	}

	## Plot the genomic distribution.
	p <- ggplot(genomic_dist, aes(x = samples, y = count, fill = fct_rev(simple_annotations))) +
		geom_col(position = "fill") +
		scale_fill_viridis_d(direction = -1, name="Annotation") +
		coord_flip() +
		ylab("Fraction") +
		theme_bw() +
		theme(
			axis.title.y = element_blank(),
			panel.grid = element_blank()
		)

	if (!is.na(metadata(genomic_distribution)$quantiles)) p <- p + facet_grid(fct_rev(factor(ntile)) ~ .)

	return(p)
}
