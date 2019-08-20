
#' Dinucleotide Analysis
#'
#' Analysis of -1 and +1 dinucleotide frequencies.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges resize score
#' @importFrom dplyr filter count mutate bind_rows group_by ungroup
#' @importFrom Rsamtools FaFile getSeq
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param samples Either 'all' or vector of sample names to analyze
#' @param genome_assembly fasta file of genome assembly
#' @param threshold TSS read threshold
#'
#' @return tibble with dinucleotide frequencies
#' @rdname dinucleotide_frequencies-function

dinucleotide_frequencies <- function(experiment, genome_assembly, samples = "all", threshold = 1) {
	## Loading genome assembly.
	fasta_assembly <- FaFile(genome_assembly)

	## Filter TSSs and get sequences.
	if (samples == "all") samples <- names(experiment@experiment$TSSs)

	sequences <- experiment@experiment$TSSs[samples] %>%
		map(~ .[score(.) >= threshold] %>%
			resize(., width=2, fix="end") %>%
			getSeq(fasta_assembly, .)
		)

	## Get dinucleotide frequencies.
	freqs <- sequences %>%
		map(~ as.character(.) %>% enframe(name = "chr", value = "dinucleotide")) %>%
		bind_rows(.id = "sample") %>%
		group_by(sample) %>%
		count(dinucleotide, sort=TRUE, name="occurences") %>%
		mutate("frequency"=occurences/sum(occurences)) %>%
		ungroup

	return(freqs)
}

#' Plot Dinucleotide Frequencies
#'
#' Plot results from dinucleotide analysis
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate pull group_by summarize arrange desc
#' @importFrom forcats fct_rev fct_reorder
#' @importFrom magrittr %>%
#'
#' @param dinucleotide_frequencies tibble from dinucleotide_frequencies analysis
#' @param ncol Number of columns to plot
#'
#' @return ggplot2 object of dinucleotide frequencies plot
#'
#' @export
#' @rdname plot_dinucleotide_frequencies-function

plot_dinucleotide_frequencies <- function(dinucleotide_frequencies, ncol = 1) {
	## Set factor order for dinucleotide.
	factor_order <- dinucleotide_frequencies %>%
		group_by(dinucleotide) %>%
		summarize(avg_freq = mean(frequency)) %>%
		arrange(desc(avg_freq)) %>%
		mutate(dinucleotide = fct_reorder(dinucleotide, desc(avg_freq))) %>%
		pull(dinucleotide) %>%
		levels
		

	dinucleotide_frequencies <- mutate(
		dinucleotide_frequencies,
		"dinucleotide" = factor(dinucleotide, levels = factor_order)
	)

	## Plot dinucleotide frequencies.
	p <- ggplot(dinucleotide_frequencies, aes(x = fct_rev(dinucleotide), y = frequency)) +
		geom_col(width = 0.5, aes(fill = frequency)) +
		theme_bw() +
		scale_fill_viridis_c(name = "Frequency") +
		coord_flip() +
		facet_wrap(~ sample, ncol = ncol) +
		labs(
			x="Dinucleotide",
			y="Frequency"
		)

	return(p)
}
