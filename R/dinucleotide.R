
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
#' @param dominant Consider only dominant
#' @param data_conditions Condition the data (filtering and quantile/grouping available)
#'
#' @return tibble with dinucleotide frequencies
#'
#' @rdname dinucleotide_frequencies-function
#'
#' @export

dinucleotide_frequencies <- function(
	experiment, genome_assembly, samples = "all",
	threshold = 1, dominant = FALSE, data_conditions = NA
) {

	## Loading genome assembly.
	fasta_assembly <- FaFile(genome_assembly)

	## Getting appropriate samples.
	select_samples <- extract_counts(experiment, "tss", samples)

	## Preliminary filtering of data.
	select_samples <- preliminary_filter(select_samples, dominant, threshold)

	## Apply conditions to data.
	if (!is.na(data_conditions)) {
		select_samples <- do.call(group_data, c(list(signal_data = select_samples), data_conditions))
	}

	## Prepare samples for analysis.
	select_samples <- rbindlist(select_samples, idcol = "sample")
	select_samples[, tss := start]
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

	## Find dinucleotide frequencies.
	groupings <- any(names(data_conditions) %in% c("quantile_by", "grouping"))

	if (!groupings) {
		freqs <- seqs[,
			.(count = .N), by = .(sample, dinucleotide)
		][,
			.(dinucleotide, count, freqs = count / sum(count)),
			by = sample
		]
	} else {
		freqs <- seqs[,
			.(count = .N), by = .(sample, dinucleotide, grouping)
		][,
			.(dinucleotide, count, freqs = count / sum(count)),
			by = .(sample, grouping)
		]
	}

	## Prepare DataFrame to return.
	freqs_df <- DataFrame(freqs)
	metadata(freqs_df)$groupings <- groupings
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
	groupings <- metadata(dinucleotide_frequencies)$groupings

	## Convert DataFrame to data.table.
	freqs <- as.data.table(dinucleotide_frequencies)
	
	## Set factor order for dinucleotide.
	if (!groupings) {
		freqs <- freqs[,
			.(sample, count, freqs, mean_freqs = mean(freqs)),
			by = dinucleotide
		]
	} else {
		freqs <- freqs[,
			.(sample, count, freqs, grouping, mean_freqs = mean(freqs)),
			by = dinucleotide
		]
	}

	freqs[,
		rank := dense_rank(mean_freqs)
	][,
		dinucleotide := fct_reorder(factor(dinucleotide), rank)
	][,
		c("mean_freqs", "rank") := NULL
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

	if (groupings) {
		p <- p + facet_wrap(grouping ~ sample)
	} else {
		p <- p + facet_wrap(~ sample, ncol = ncol)
	}

	return(p)
}
