
#' Genomic Distribution
#'
#' Get genomic distribution of TSSs and TSRs
#'
#' @import tibble
#' @importFrom dplyr mutate count case_when arrange rename
#'
#' @param experiment tsrexplorer object with annotated TSRs
#' @param sample Name of sample to analyze
#' @param data_type Whether to get distribution of TSSs or TSRs
#' @param threshold Filter out TSSs or TSRs under a certain read count number
#'
#' @return tibble with TSS or TSR genomic distribution stats
#'
#' @export 
#' @rdname genomic_distribution-function

genomic_distribution <- function(experiment, sample, data_type = c("tss", "tsr"), threshold = 1) {

	## Pull data from tsrexplorer object.
	if (data_type == "tss") {
		selected_sample <- experiment@annotated$TSSs[[sample]]
	} else if (data_type == "tsr") {
		selected_sample <- experiment@annotated$TSRs[[sample]] %>%
			rename("score" = nTAGs)
	}

	## Prepare data
	genomic_distribution <- selected_sample %>%
		filter(score >= threshold) %>%
		mutate(annotation = case_when(
			annotation == "Promoter" ~ "Promoter",
			grepl(annotation, pattern="Exon") ~ "Exon",
			grepl(annotation, pattern="Intron") ~ "Intron",
			grepl(annotation, pattern="Downstream") ~ "Downstream",
			annotation == "Distal Intergenic" ~ "Intergenic"
		)) %>%
		count(annotation, name = "count") %>%
		mutate(
			annotation = factor(annotation, levels = c(
				"Promoter", "Exon", "Intron", "Downstream", "Intergenic"
			)),
			fraction = count / sum(count)
		) %>%
		arrange(annotation)

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
#' @export 
#' @rdname plot_genomic_distribution-function

plot_genomic_distribution <- function(genomic_distribution) {
	
	## Prepare data for plotting.
	genomic_distribution <- mutate(genomic_distribution, sample = "samp1")
	
	## Plot genomic distribution.
	p <- ggplot(genomic_distribution, aes(x=sample, y=count, fill=fct_rev(annotation))) +
		geom_col(position="fill") +
		scale_fill_viridis_d(direction=-1, name="Annotation") +
		coord_flip() +
		ylab("Fraction") +
		theme_minimal() +
		theme(
			axis.text.y=element_blank(),
			axis.title.y=element_blank(),
			panel.grid=element_blank()
		)
}
