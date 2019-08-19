
#' Genomic Distribution
#'
#' Get genomic distribution of TSSs and TSRs
#'
#' @import tibble
#' @importFrom dplyr mutate count case_when arrange rename bind_rows group_by ungroup
#'
#' @param experiment tsrexplorer object with annotated TSRs
#' @param samples Either "all" or vector of sample names
#' @param data_type Whether to get distribution of TSSs or TSRs
#' @param threshold Filter out TSSs or TSRs under a certain read count number
#'
#' @return tibble with TSS or TSR genomic distribution stats
#'
#' @export 
#' @rdname genomic_distribution-function

genomic_distribution <- function(experiment, data_type = c("tss", "tsr"), samples = "all", threshold = 1) {

	## Pull data from tsrexplorer object.
	if (data_type == "tss") {
		if (samples == "all") samples <- names(experiment@annotated$TSSs)
		selected_samples <- experiment@annotated$TSSs[samples] %>%
			bind_rows(.id = "samples")
	} else if (data_type == "tsr") {
		if (samples == "all") samples <- names(experiment@annotated$TSRs)
		selected_samples <- experiment@annotated$TSRs[samples] %>%
			bind_rows(.id = "samples") %>%
			rename("score" = nTAGs)
	}

	## Prepare data
	genomic_distribution <- selected_samples%>%
		filter(score >= threshold) %>%
		mutate(annotation = case_when(
			annotation == "Promoter" ~ "Promoter",
			grepl(annotation, pattern="Exon") ~ "Exon",
			grepl(annotation, pattern="Intron") ~ "Intron",
			grepl(annotation, pattern="Downstream") ~ "Downstream",
			annotation == "Distal Intergenic" ~ "Intergenic"
		)) %>%
		count(samples, annotation, name = "count") %>%
		mutate(
			annotation = factor(annotation, levels = c(
				"Promoter", "Exon", "Intron", "Downstream", "Intergenic"
		))) %>%
		group_by(samples) %>%
		mutate(fraction = count / sum(count)) %>%
		ungroup %>%
		arrange(samples, annotation)

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
#' @param ncol Number of columns when plotting
#'
#' @return ggplot2 object with TSS or TSR genomic distribution plot
#'
#' @export 
#' @rdname plot_genomic_distribution-function

plot_genomic_distribution <- function(genomic_distribution, ncol = 1) {
	
	## Plot genomic distribution.
	p <- ggplot(genomic_distribution, aes(x=samples, y=count, fill=fct_rev(annotation))) +
		geom_col(position="fill") +
		scale_fill_viridis_d(direction=-1, name="Annotation") +
		coord_flip() +
		ylab("Fraction") +
		theme_minimal() +
		theme(
			axis.title.y = element_blank(),
			panel.grid = element_blank()
		)
}
