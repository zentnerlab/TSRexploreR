
#' Genomic Distribution
#'
#' Get genomic distribution of TSSs and TSRs
#'
#' @import tibble
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
#'
#' @return tibble with TSS or TSR genomic distribution stats
#'
#' @rdname genomic_distribution-function
#'
#' @export 

genomic_distribution <- function(experiment, data_type = c("tss", "tsr"), samples = "all", threshold = 1, quantiles = NA) {

	## Pull data from tsrexplorer object.
	if (data_type == "tss") {
		if (samples == "all") samples <- names(experiment@counts$TSSs$raw)
		selected_samples <- extract(experiment@counts$TSSs$raw, samples)
	} else if (data_type == "tsr") {
		if (samples == "all") samples <- names(experiment@counts$TSRs$raw)
		selected_samples <- extract(experiment@counts$TSRs$raw, samples)
	}


	## Initial preparation of data.
	selected_samples <- selected_samples %>%
		map(function(x) {
			ranges <- rowRanges(x) %>% as_tibble(.name_repair = "unique")
			scores <- assay(x, "raw") %>% as.numeric
			
			ranges <- mutate(ranges, score = scores)
			return(ranges)
		}) %>%
		bind_rows(.id = "samples") %>%
		filter(score >= threshold)

	## Break data into quantiles if quantiles set is greater than 1.
	if (!is.na(quantiles)) {
		selected_samples <- selected_samples %>%
			group_by(samples) %>%
			mutate(ntile = ntile(score, quantiles))
	}
	
	## Clean up genomic annotations.
	genomic_distribution <- selected_samples %>%
		mutate(annotation = case_when(
			annotation == "Promoter" ~ "Promoter",
			grepl(annotation, pattern="(Exon|UTR)") ~ "Exon",
			grepl(annotation, pattern="Intron") ~ "Intron",
			grepl(annotation, pattern="Downstream") ~ "Downstream",
			annotation == "Distal Intergenic" ~ "Intergenic"
		)) %>%
		mutate(
			annotation = factor(annotation, levels = c(
				"Promoter", "Exon", "Intron", "Downstream", "Intergenic"
			))
		)

	## Prepade data to be plotted later.
	if (!is.na(quantiles)) {
		genomic_distribution <- genomic_distribution %>%
			count(samples, annotation, ntile, name = "count") %>%
			group_by(samples, ntile) %>%
			mutate(fraction = count / sum(count)) %>%
			ungroup %>%
			arrange(samples, ntile, annotation)
	} else {
		genomic_distribution <- genomic_distribution %>%
			count(samples, annotation, name = "count") %>%
			group_by(samples) %>%
			mutate(fraction = count / sum(count)) %>%
			ungroup %>%
			arrange(samples, annotation)
	}

	## Prepare summarized experiment to return.
	genome_fracs <- genomic_distribution %>%
		select(-count) %>%
		spread(samples, fraction)
	if (!is.na(quantiles)) {
		genome_fracs <- genome_fracs %>%
			select(-annotation, -ntile) %>%
			as.matrix
	} else {
		genome_fracs <- genome_fracs %>%
			select(-annotation) %>%
			as.matrix
	}

	genome_counts <- genomic_distribution %>%
		select(-fraction) %>%
		spread(samples, count)
	if (!is.na(quantiles)) {
		genome_counts <- genome_counts %>%
			select(-annotation, -ntile) %>%
			as.matrix
	} else {
		genome_counts <- genome_counts %>%
			select(-annotation, -ntile) %>%
			as.matrix
	}

	if (!is.na(quantiles)) {
		row_data <- distinct(genomic_distribution, annotation, ntile)
	} else {
		row_data <- distinct(genomic_distribution, annotation)
	}

	col_data <- DataFrame(sample = colnames(genome_fracs))

	dist_exp <- SummarizedExperiment(
		assay = list(counts = genome_counts, fractions = genome_fracs),
		rowData = row_data,
		colData = col_data
	)

	## Add quantile information to summarized experiment.
	if (!is.na(quantiles)) {
		metadata(dist_exp)$ntiles <- TRUE
	} else {
		metadata(dist_exp)$ntiles <- FALSE
	}

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
	
	## Format summarized experiment for plotting.
	counts <- assay(genomic_distribution, "counts") %>%
		as_tibble(.name_repair = "unique")
	annotations <- rowData(genomic_distribution) %>%
		as_tibble(.name_repair = "unique")

	genomic_dist <- bind_cols(annotations, counts) %>%
		gather(-annotation, key = "samples", value = "counts")

	## Order the samples if required.
	if (!is.na(sample_order)) {
		genomic_dist <- genomic_dist %>%
			mutate(samples = fct_rev(factor(samples, levels = sample_order)))
	}

	## Plot the genomic distribution.
	p <- ggplot(genomic_dist, aes(x = samples, y = counts, fill = fct_rev(annotation))) +
		geom_col(position = "fill") +
		scale_fill_viridis_d(direction = -1, name="Annotation") +
		coord_flip() +
		ylab("Fraction") +
		theme_bw() +
		theme(
			axis.title.y = element_blank(),
			panel.grid = element_blank()
		)

	if (metadata(genomic_distribution)$ntiles) p <- p + facet_grid(fct_rev(factor(ntile)) ~ .)

	return(p)
}
