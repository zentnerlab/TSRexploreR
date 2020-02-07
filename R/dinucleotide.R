
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

dinucleotide_frequencies <- function(experiment, genome_assembly, samples = "all", threshold = 1, quantiles = NA) {

	## Loading genome assembly.
	fasta_assembly <- FaFile(genome_assembly)

	## Getting appropriate samples.
	select_samples <- experiment %>%
		extract_counts("tss", samples) %>%
		bind_rows(.id = "sample") %>%
		as.data.table

	## Prepare samples for analysis.
	select_samples <- select_samples[
		score >= threshold,
		.(sample, seqnames, start, end, strand, score, tss = start)
	]
	select_samples[,
		c("start", "end") := list(
			ifelse(strand == "+", start - 1, start),
			ifelse(strand == "+", end, end + 1)
		)
	]

	## Get dinucleotides.
	seqs <- select_samples %>%
		makeGRangesFromDataFrame %>%
		getSeq(fasta_assembly, .) %>%
		as.data.table %>%
		bind_cols(select_samples, .)
	setnames(seqs, old = "x", new = "dinucleotide")

	## Add ntile information if required.
	if (!is.na(quantiles)) {
		seqs[, ntile := ntile(score, quantiles), by = sample]
	}

	## Find dinucleotide frequencies.
	if (is.na(quantiles)) {
		freqs <- seqs[,
			.(count = .N), by = .(sample, dinucleotide)
		][,
			.(dinucleotide, count, freqs = count / sum(count)),
			by = sample
		]
	} else {
		freqs <- seqs[,
			.(count = .N), by = .(sample, dinucleotide, ntile)
		][,
			.(dinucleotide, count, freqs = count / sum(count)),
			by = .(sample, ntile)
		]
	}

	## Prepare DataFrame to return.
	freqs_df <- DataFrame(freqs)
	metadata(freqs_df)$quantiles <- quantiles
	metadata(freqs_df)$threshold <- threshold
	
	return(freqs_df)
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
#' @param sample_order Optional vector of sample names to order plots
#' @param ... Arguments passed to geom_col
#'
#' @return ggplot2 object of dinucleotide frequencies plot
#'
#' @rdname plot_dinucleotide_frequencies-function
#'
#' @export

plot_dinucleotide_frequencies <- function(dinucleotide_frequencies, sample_order = NULL, ncol = 1, ...) {
	
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

	## Change sample order if specified.
	if (!is.null(sample_order)) {
		frequencies <- mutate(frequencies, "sample" = factor(sample, levels = sample_order))
	}

	## Plot dinucleotide frequencies.
	p <- ggplot(frequencies, aes(x = fct_rev(dinucleotide), y = frequency)) +
		geom_col(width = 0.5, aes(fill = frequency), ...) +
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
