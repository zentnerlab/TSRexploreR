
#' Dinucleotide Analysis
#'
#' Analysis of -1 and +1 dinucleotide frequencies.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges resize score
#' @importFrom dplyr mutate bind_rows group_by ntile ungroup select count
#' @importFrom Rsamtools FaFile getSeq
#' @importFrom purrr map map2
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param samples Either 'all' or vector of sample names to analyze
#' @param genome_assembly fasta file of genome assembly
#' @param threshold TSS read threshold
#' @param quantiles Number of quantiles to break data into
#'
#' @return tibble with dinucleotide frequencies
#'
#' @rdname dinucleotide_frequencies-function
#'
#' @export

dinucleotide_frequencies <- function(experiment, genome_assembly, samples = "all", threshold = 1, quantiles = 1) {
	## Loading genome assembly.
	fasta_assembly <- FaFile(genome_assembly)

	## Filter TSSs and get sequences.
	if (samples == "all") samples <- names(experiment@counts$TSSs$raw)
	select_samples <- experiment@counts$TSSs$raw[samples]
	
	## Prepare samples for analysis.
	select_samples

	select_samples <- map(select_samples, ~ .[score(.) >= threshold])

	sample_scores <- map(select_samples, ~ score(.))

	sequences <- select_samples %>%
		map(
			~ resize(., width=2, fix="end") %>%
				getSeq(fasta_assembly, .)
		)
			

	## Get dinucleotide frequencies.
	freqs <- sequences %>%
		map(~ as.character(.) %>% enframe(name = NULL, value = "dinucleotide")) %>%
		map2(., sample_scores, ~ mutate(.x, score = .y)) %>%
		bind_rows(.id = "sample")

	if (quantiles > 1) {
		freqs <- freqs %>%
			group_by(sample) %>%
			mutate(ntile = ntile(score, quantiles)) %>%
			ungroup %>%
			select(-score) %>%
			count(sample, dinucleotide, ntile, name = "occurences") %>%
			group_by(sample, ntile) %>%
			mutate("frequency" = occurences/sum(occurences)) %>%
			ungroup

		freqs <- list(
			"quantile_plot" = TRUE,
			"dinucleotide_frequencies" = freqs
		)
	} else {
		freqs <- freqs %>%
			count(sample, dinucleotide, name = "occurences") %>%
			group_by(sample) %>%
			mutate("frequency" = occurences/sum(occurences)) %>%
			ungroup

		freqs <- list(
			"quantile_plot" = FALSE,
			"dinucleotide_frequencies" = freqs
		)
	}

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
#' @param ncol Number of columns to plot if not quantile plot
#'
#' @return ggplot2 object of dinucleotide frequencies plot
#'
#' @rdname plot_dinucleotide_frequencies-function
#'
#' @export

plot_dinucleotide_frequencies <- function(dinucleotide_frequencies, ncol = 1) {
	
	## Set factor order for dinucleotide.
	factor_order <- dinucleotide_frequencies$dinucleotide_frequencies %>%
		group_by(dinucleotide) %>%
		summarize(avg_freq = mean(frequency)) %>%
		arrange(desc(avg_freq)) %>%
		mutate(dinucleotide = fct_reorder(dinucleotide, desc(avg_freq))) %>%
		pull(dinucleotide) %>%
		levels
		

	frequencies <- mutate(
		dinucleotide_frequencies$dinucleotide_frequencies,
		"dinucleotide" = factor(dinucleotide, levels = factor_order)
	)

	## Plot dinucleotide frequencies.
	p <- ggplot(frequencies, aes(x = fct_rev(dinucleotide), y = frequency)) +
		geom_col(width = 0.5, aes(fill = frequency)) +
		theme_bw() +
		scale_fill_viridis_c(name = "Frequency") +
		coord_flip() +
		labs(
			x="Dinucleotide",
			y="Frequency"
		)

	if (dinucleotide_frequencies$quantile_plot) {
		p <- p + facet_grid(fct_rev(factor(ntile)) ~ sample)
	} else {
		p <- p + facet_wrap(~ sample, ncol = ncol)
	}

	return(p)
}
