
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
#' @param dominant Consider only dominant
#'
#' @return tibble with dinucleotide frequencies
#'
#' @rdname dinucleotide_frequencies-function
#'
#' @export

dinucleotide_frequencies <- function(
	experiment, genome_assembly, samples = "all",
	threshold = 1, quantiles = NA, dominant = FALSE
) {

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
		.(sample, seqnames, start, end, strand, dominant,  score, tss = start)
	]
	select_samples[,
		c("start", "end") := list(
			ifelse(strand == "+", start - 1, start),
			ifelse(strand == "+", end, end + 1)
		)
	]

	## Consider only dominant TSSs if required.
	if (dominant) {
		select_samples <- select_samples[(dominant)]
	}
	select_samples[, dominant := NULL]

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

plot_dinucleotide_frequencies <- function(
	dinucleotide_frequencies, sample_order = NA, ncol = 1, ...) {

	## Pull out some info from the DataFrame.
	quantiles <- metadata(dinucleotide_frequencies)$quantiles

	## Convert DataFrame to data.table.
	freqs <- as.data.table(dinucleotide_frequencies)
	
	## Set factor order for dinucleotide.
	if (is.na(quantiles)) {
		freqs <- freqs[,
			.(sample, count, freqs, mean_freqs = mean(freqs)),
			by = dinucleotide
		]
	} else {
		freqs <- freqs[,
			.(sample, count, freqs, ntile, mean_freqs = mean(freqs)),
			by = dinucleotide
		]
	}

	freqs[,
		rank := dense_rank(mean_freqs)
	][,
		dinucleotide := fct_reorder(factor(dinucleotide), rank)
	][,
		c("mean_freqs", "rank") := list(NULL, NULL)
	]

	## Change sample order if specified.
	if (!is.na(sample_order)) {
		freqs[, sample := factor(sample, levels = sample_order)]
	}

	## Plot dinucleotide frequencies.
	p <- ggplot(freqs, aes(x = dinucleotide, y = freqs)) +
		geom_col(width = 0.5, aes(fill = freqs), ...) +
		theme_bw() +
		scale_fill_viridis_c(name = "Frequency") +
		coord_flip() +
		labs(
			x="Dinucleotide",
			y="Frequency"
		)

	if (!is.na(quantiles)) {
		p <- p + facet_grid(fct_rev(factor(ntile)) ~ sample)
	} else {
		p <- p + facet_wrap(~ sample, ncol = ncol)
	}

	return(p)
}
