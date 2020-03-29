
#' Retrieve Sequences Near TSSs
#'
#' Retrieve sequences surrounding TSSs for further plotting
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @import data.table
#' @importFrom S4Vectors DataFrame "metadata<-"
#' @importFrom Rsamtools FaFile
#' @importFrom purrr map map2
#' @importFrom dplyr mutate filter ntile pull bind_rows
#' @importFrom Biostrings getSeq DNAStringSet
#' @importFrom GenomicRanges GRanges resize width makeGRangesFromDataFrame
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param samples Either "all" or names of samples to analyze
#' @param genome_assembly Genome assembly in fasta format or bioconductor BSgenome
#' @param threshold Keep only TSSs with threshold number of reads or more
#' @param distance Bases to add on each side of TSS
#' @param dominant Whether only dominant should be considered
#' @param data_condtions Condition the data (filter, quantile, and grouping available)
#' @return Sequences surrounding TSSs
#'
#' @rdname tss_sequences-function
#'
#' @export

tss_sequences <- function(
	experiment, samples = "all", genome_assembly, threshold = 1,
	distance = 10, dominant = FALSE, data_conditions = NA
) {

	## Open genome assembly.
	if (is(genome_assembly, "character")) {
		genome_assembly <- FaFile(genome_assembly)
	} else if (is(genome_assembly, "BSgenome")) {
		genome_assembly <- genome_assembly
	}

	## Pull selected samples.
	select_samples <- extract_counts(experiment, "tss", samples)

	## Preliminary filtering of data.
	select_samples <- preliminary_filter(select_samples, dominant, threshold)

	## Add conditions to data.
	if (!is.na(data_conditions)) {
		select_samples <- do.call(group_data, c(list(signal_data = select_samples), data_conditions))
	}

	## Expand ranges.
	select_samples <- rbindlist(select_samples, id = "sample")
	select_samples[, c("start", "end", "tss") := list(start - distance, end + distance, start)]

	## Get sequences.
	seqs <- select_samples %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
		getSeq(genome_assembly, .) %>%
		as.data.table %>%
		bind_cols(select_samples, .)
			
	setnames(seqs, old = "x", new = "sequence")
	
	## Generate return DataFrame
	groupings <- any(names(data_conditions) %in% c("quantile_by", "grouping"))
	
	seqs <- DataFrame(seqs)
	metadata(seqs)$groupings <- groupings
	metadata(seqs)$threshold <- threshold
	metadata(seqs)$distance <- distance
	metadata(seqs)$dominant <- dominant

	return(seqs)
}

#' Generate Sequence Logo
#'
#' Create a sequence logo for the sequences around TSSs
#'
#' @import tibble
#' @import ggplot2
#' @import ggseqlogo
#' @importFrom Biostrings DNAStringSet consensusMatrix
#' @importFrom purrr map pmap
#' @importFrom dplyr pull
#' @importFrom tidyr unite separate
#' @importFrom magrittr %>%
#' @importFrom cowplot plot_grid
#'
#' @param tss_sequences Sequences surrounding TSS generated with tss_sequences
#' @param ncol Number of columns to plot if quantiles is not set
#' @param font_size Font size plots
#'
#' @return ggplot2 object with sequence logo
#'
#' @rdname plot_sequence_logo-function
#'
#' @export

plot_sequence_logo <- function(tss_sequences, ncol = 1, font_size = 6) {

	## Get some info used to pull sequencs.
	distance <- metadata(tss_sequences)$distance
	groupings <- metadata(tss_sequences)$groupings

	## Grab sequences from input.
	if (!groupings) {
		sequences <- tss_sequences %>%
			as.data.table %>%
			split(.$sample) %>%
			map(function(x) {x[["sequence"]]})
	} else {
		sequences <- tss_sequences %>%
			as.data.table %>%
			split(.$grouping) %>%
			map(function(x) {
				split(x, x$sample) %>%
				map(function(y) {y[["sequence"]]})
			})
	}

	## Create viridis color scheme for bases.
	viridis_bases <- make_col_scheme(
		chars = c("A", "C", "G", "T"),
		groups = c("A", "C", "G", "T"),
		cols = c("#431352", "#34698c", "#44b57b", "#fde540")
	)

	## Make sequence logo.
	if (!groupings) {
		p <- ggseqlogo(sequences, ncol = ncol) +
			theme(text = element_text(size = font_size)) #+
			#scale_x_continuous(
			#	breaks = c(1, distance, distance + 1, (distance * 2) + 1),
			#	labels = c(-distance, -1, +1, distance + 1)
			#)
	} else {
		p <- rev(sequences) %>%
			map(function(x) {
				ggseqlogo(x, ncol = ncol) +
					theme(text = element_text(size = font_size)) #+
					#scale_x_continuous(
					#	breaks = c(1, distance, distance + 1, (distance * 2) + 1),
					#	labels = c(-distance, -1, +1, distance + 1)
					#)
			})

		p <- plot_grid(plotlist = p, labels = rev(names(sequences)), ncol = 1)
	}

	return(p)
}

#' Plot Sequence Colormap
#'
#' Make a sequence colormap for the sequences around TSSs.
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr bind_rows rename pull bind_cols mutate group_by ungroup desc
#' @importFrom purrr map
#' @importFrom stringr str_split
#' @importFrom tidyr separate gather
#' @importFrom forcats fct_reorder fct_rev
#' @importFrom Biostrings DNAStringSet
#' @importFrom BiocGenerics width
#'
#' @param tss_sequences Sequences surrounding TSS generated with tss_sequences
#' @param ncol Number of columns to plot data if quantiles not specified
#' @param base_colors Named vector specifying colors for each base
#' @param text_size Size of text for plots
#' @param ... Arguments passed to geom_tile
#'
#' @return ggplot2 object of sequence colormap
#'
#' @rdname plot_sequence_colormap-function
#'
#' @export

plot_sequence_colormap <- function(
	tss_sequences,
	ncol = 1,
	base_colors = c(
		"A" = "#109649",
		"C" = "#255C99",
		"G" = "#F7B32C",
		"T" = "#D62839"
	),
	text_size = 6,
	...
) {
	## Grab some information out of DataFrame.
	distance <- metadata(tss_sequences)$distance
	quantiles <- metadata(tss_sequences)$quantiles

	## Start preparing data for plotting.
	seq_data <- as.data.table(tss_sequences)
	seq_data[, width := nchar(sequence)]

	## Split sequences into columns
	split_seqs <- seq_data[, as.data.table(str_split(sequence, "", simplify = TRUE))]
	setnames(
		split_seqs,
		old = seq(1, (distance * 2) + 1),
		new = as.character(c(seq(-distance, -1), seq(1, distance + 1)))
	)
	seq_data <- bind_cols(seq_data, split_seqs)

	## Get order of TSSs for plotting.
	if (!is.na(quantiles)) {
		seq_data[, rank := dense_rank(score), by = .(sample, ntile)]
		seq_data[, name := sprintf("TSS%010d", seq_len(.N)), by = .(sample, ntile)]
	} else {
		seq_data[, rank := dense_rank(score), by = sample]
		seq_data[, name := sprintf("TSS%010d", seq_len(.N)), by = sample]
	}

	## Format data for plotting.
	long_data <- seq_data %>%
		melt(
			measure.vars = as.character(c(seq(-distance, -1), seq(1, distance + 1))),
			variable.name = "position", value.name = "base"
		)

	long_data[,
		c("name", "position", "base") := list(
			name = fct_reorder(name, rank),
			position = as.numeric(position),
			base = factor(base, levels = c("A", "C", "G", "T"))
		)
	]
		
	## Plot sequence colormap
	p <- ggplot(long_data, aes(x = position, y = name)) +
		geom_tile(aes(fill = base, color = base)) +
		scale_fill_manual(values = base_colors) +
		scale_color_manual(values = base_colors) +
		theme_minimal() +
		theme(
			axis.title.y = element_blank(),
			axis.text.y = element_blank(),
			legend.title = element_blank(),
			axis.title.x = element_text(margin = margin(t = 15)),
			panel.grid = element_blank(),
			text = element_text(size = text_size)
		) +
		scale_x_continuous(
			breaks = c(1, distance, distance + 1, (distance * 2) + 1),
			labels = c(-distance, -1, 1, distance + 1)
		)

	if (is.na(quantiles)) {
		p <- p + facet_wrap(. ~ sample, scales = "free", ncol = ncol)
	} else {
		p <- p + facet_wrap(ntile ~ sample, scales = "free", ncol = ncol)
	}

	return(p)
}
