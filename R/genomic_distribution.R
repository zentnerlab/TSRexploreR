
#' TSR Genomic Distribution
#'
#' Get genomic distribution of TSRs
#'
#' @import tibble
#' @importFrom dplyr mutate count case_when
#'
#' @param experiment tsrexplorer object with annotated TSRs
#' @param sample Name of sample to analyze
#'
#' @return tibble with TSR genomic distribution stats
#'
#' @export 
#' @rdname tsr_genomic_distribution-function

tsr_genomic_distribution <- function(experiment, sample) {
	## Prepare data
	genomic_distribution <- experiment@annotated$TSRs[[sample]] %>%
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

#' Plot TSR Genomic Distribution
#'
#' Get genomic distribution of TSRs
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom forcats fct_rev
#'
#' @param tsr_genomic_distribution tibble of TSR genomic distributions from tsr_genomic_distribution
#'
#' @return ggplot2 object with TSR genomic distribution plot
#'
#' @export 
#' @rdname plot_tsr_genomic_distribution-function

plot_tsr_genomic_distribution <- function(tsr_genomic_distribution) {
	## Prepare data for plotting
	tsr_genomic_distribution <- mutate(tsr_genomic_distribution, sample = "samp1")
	
	## Plot genomic distribution
	p <- ggplot(tsr_genomic_distribution, aes(x=sample, y=count, fill=fct_rev(annotation))) +
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
