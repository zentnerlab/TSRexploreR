
#' Dinucleotide Analysis
#'
#' Analysis of -1 and +1 dinucleotide frequencies.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges resize score
#' @importFrom dplyr filter count mutate
#' @importFrom Rsamtools FaFile getSeq
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param sample Name of sample to analyze
#' @param genome_assembly fasta file of genome assembly
#' @param threshold TSS read threshold
#'
#' @return tibble with dinucleotide frequencies
#' @rdname dinucleotide_frequencies-function

dinucleotide_frequencies <- function(experiment, sample, genome_assembly, threshold=1) {
	## Loading genome assembly.
	assembly <- FaFile(genome_assembly)

	## Filter TSSs and get sequences.
	sequences <- experiment@experiment$TSSs[[sample]] %>%
		.[score(.) >= threshold] %>%
		resize(., width=2, fix="end") %>%
		getSeq(assembly, .)

	## Get dinucleotide frequencies.
	freqs <- sequences %>%
		as.character %>%
		enframe(name="chr", value="dinucleotide") %>%
		count(dinucleotide, sort=TRUE, name="occurences") %>%
		mutate("frequency"=occurences/sum(occurences))

	return(freqs)
}

#' Plot Dinucleotide Frequencies
#'
#' Plot results from dinucleotide analysis
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate pull
#' @importFrom forcats fct_rev
#' @importFrom magrittr %>%
#'
#' @param dinucleotide_frequencies tibble from dinucleotide_frequencies analysis
#'
#' @return ggplot2 object of dinucleotide frequencies plot
#'
#' @export
#' @rdname plot_dinucleotide_frequencies-function

plot_dinucleotide_frequencies <- function(dinucleotide_frequencies) {
	## Set factor order for dinucleotide.
	dinucleotide_frequencies <- mutate(
		dinucleotide_frequencies,
		"dinucleotide"=factor(dinucleotide, levels=pull(dinucleotide_frequencies, dinucleotide) %>% unique)
	)

	## Plot dinucleotide frequencies.
	p <- ggplot(dinucleotide_frequencies, aes(x=fct_rev(dinucleotide), y=frequency)) +
		geom_col(width=0.5, aes(fill=frequency)) +
		theme_classic() +
		scale_fill_viridis_c(name="Frequency") +
		coord_flip() +
		labs(
			x="Dinucleotide",
			y="Frequency"
		)

	return(p)
}
