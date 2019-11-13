
#' Genomic Distribution
#'
#' Get genomic distribution of TSSs and TSRs
#'
#' @import tibble
#' @importFrom dplyr mutate count case_when arrange rename bind_rows group_by ungroup ntile
#'
#' @param experiment tsrexplorer object with annotated TSRs
#' @param samples Either "all" or vector of sample names
#' @param data_type Whether to get distribution of TSSs or TSRs
#' @param threshold Filter out TSSs or TSRs under a certain read count number
#' @param number of quantiles to break the data down into
#'
#' @return tibble with TSS or TSR genomic distribution stats
#'
#' @rdname genomic_distribution-function
#'
#' @export 

genomic_distribution <- function(experiment, data_type = c("tss", "tsr"), samples = "all", threshold = 1, quantiles = 1) {

	## Pull data from tsrexplorer object.
	if (data_type == "tss") {
		if (samples == "all") samples <- names(experiment@annotated$TSSs$raw)
		selected_samples <- experiment@annotated$TSSs$raw[samples]
	} else if (data_type == "tsr") {
		if (samples == "all") samples <- names(experiment@annotated$TSRs$raw)
		selected_samples <- experiment@annotated$TSRs$raw[samples]
	}


	## Initial preparation of data.
	selected_samples <- selected_samples %>%
		bind_rows(.id = "samples") %>%
		filter(score >= threshold)

	## Break data into quantiles if quantiles set is greater than 1.
	if (quantiles > 1) {
		selected_samples <- selected_samples %>%
			group_by(samples) %>%
			mutate(ntile = ntile(score, quantiles))
	}
	
	## Clean up genomic annotations.
	genomic_distribution <- selected_samples %>%
		mutate(annotation = case_when(
			annotation == "Promoter" ~ "Promoter",
			grepl(annotation, pattern="Exon") ~ "Exon",
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
	if (quantiles > 1) {
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

	## Add information on whether quantiles should be plotted
	if (quantiles > 1) {
		genomic_distribution <- list(
			"quantile_plot" = TRUE,
			"genomic_distribution" = genomic_distribution
		)
	} else {
		genomic_distribution <- list(
			"quantile_plot" = FALSE,
			"genomic_distribution" = genomic_distribution
		)
	}

	return(genomic_distribution)
}

#' Plot Genomic Distribution
#'
#' Plot genomic distribution of TSSs or TSRs
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom forcats fct_rev
#'
#' @param genomic_distribution tibble of TSS or TSR genomic distributions from tsr_genomic_distribution
#'
#' @return ggplot2 object with TSS or TSR genomic distribution plot
#'
#' @rdname plot_genomic_distribution-function
#'
#' @export 

plot_genomic_distribution <- function(genomic_distribution) {
	
	p <- ggplot(genomic_distribution$genomic_distribution, aes(x = samples, y = count, fill = fct_rev(annotation))) +
		geom_col(position = "fill") +
		scale_fill_viridis_d(direction = -1, name="Annotation") +
		coord_flip() +
		ylab("Fraction") +
		theme_bw() +
		theme(
			axis.title.y = element_blank(),
			panel.grid = element_blank()
		)

	if (genomic_distribution$quantile_plot) p <- p + facet_grid(fct_rev(factor(ntile)) ~ .)

	return(p)
}
