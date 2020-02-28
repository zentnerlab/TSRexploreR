
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
#' @param dominant Should dominant TSS per 'gene' or 'tsr' only be considered
#' @param data_group List of group options (filtering and quantiles available)
#'
#' @return tibble with TSS or TSR genomic distribution stats
#'
#' @rdname genomic_distribution-function
#'
#' @export 

genomic_distribution <- function(
	experiment, data_type = c("tss", "tsr"), samples = "all",
	threshold = NA, dominant = FALSE, data_group = NA
) {

	## Extract samples.
	selected_samples <- extract_counts(experiment, data_type, samples)

	## Preliminary sample preparation.
	if (dominant | !is.na(threshold)) {
		selected_samples <- preliminary_filter(selected_samples, dominant, threshold)
	}

	walk(selected_samples, function(x) {
        	x[, simple_annotations := factor(
                        simple_annotations,
                        levels = c("Promoter", "Exon", "Intron", "Downstream", "Intergenic")
                )]
	})


	## Apply advanced grouping.
	if (!is.na(data_group)) {
		selected_samples <- do.call(group_data, c(list(signal_data = selected_samples), data_group))
	}

	## Prepare data to be plotted later.
	selected_samples <- bind_rows(selected_samples, .id = "samples")
	quantiles <- any(names(data_group) == "quantile_by")

	if (quantiles) {
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
