
#' Retrieve Sequences Near TSSs
#'
#' Retrieve sequences surrounding TSSs for further plotting
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom Rsamtools FaFile
#' @importFrom purrr map map2
#' @importFrom dplyr mutate filter ntile pull bind_rows
#' @importFrom Biostrings getSeq DNAStringSet
#' @importFrom GenomicRanges GRanges resize width makeGRangesFromDataFrame
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param samples Either "all" or names of samples to analyze
#' @param genome_assembly Genome assembly in fasta format
#' @param threshold Keep only TSSs with threshold number of reads or more
#' @param distance Bases to add on each side of TSS
#' @param quantiles Break data into quantiles
#'
#' @return Sequences surrounding TSSs
#'
#' @rdname tss_sequences-function
#'
#' @export

tss_sequences <- function(experiment, samples = "all", genome_assembly, threshold = 1, distance = 10, quantiles = 1) {

	## Open genome assembly.
	genome_assembly <- FaFile(genome_assembly)

	## Pull selected samples.
	if (samples == "all") samples <- names(experiment@counts$TSSs$raw)
	select_samples <- experiment@counts$TSSs$raw[samples]

	## Start preparing data for retrieving sequences.
	tss_sequences <- select_samples %>%
		map(~ as_tibble(., .name_repair = "unique")) %>%
		bind_rows(.id = "sample") %>%
		filter(score >= threshold)

	tss_sequences <- mutate(tss_sequences, ntile = ntile(score, quantiles))

	## Store sequence names.
	tss_sequence_names <- tss_sequences %>%
		split(.$sample) %>%
		map(~ mutate(., tss_name = paste(seqnames, start, end, strand, score, ntile, sep = "_")))

	tss_sequence_names <- map(tss_sequence_names, ~ pull(., tss_name))

	## Get sequences for TSSs.
	tss_sequences <- tss_sequences %>%
		split(.$sample) %>%
		map(
			~ makeGRangesFromDataFrame(.) %>%
				resize(width = distance, fix = "start") %>%
				resize(width = width(.) + distance, fix = "end") %>%
				getSeq(genome_assembly, .)
		)

	## Set the names of the sequences.
	rename_sequences <- function(x, y) {
		names(x) <- y
		return(x)
	}

	tss_sequences <- map2(tss_sequences, tss_sequence_names, ~ rename_sequences(.x, .y))

	return(tss_sequences)
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
#' @importFrom dplyr mutate
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

	## Grab sequences from input.
	tss_seqs <- tss_sequences %>%
		map(~ as.data.frame(.x) %>% as_tibble(.name_repair = "unique", rownames = "position"))

	tss_seqs <- tss_seqs %>%
		map(
			~ separate(.x, position, into = c("seqnames", "start", "end", "strand", "score", "ntile"), sep = "_") %>%
				unite("position", seqnames, start, end, strand, score, sep = "_") %>%
				split(.$ntile) %>%
				map(~ pull(.x) %>% DNAStringSet)
		) %>%
		enframe

	## Calculate consensus matrix.
	consensus_matrix <- tss_seqs %>%
		pmap(function(name, value) {
			map(
				value,
				~ consensusMatrix(., as.prob = TRUE) %>%
					.[rownames(.) %in% c("A", "C", "G", "T"), ]
			) %>%
			rev
		}) %>%
		enframe %>%
		mutate("name" = pull(tss_seqs, name))

	## Create viridis color scheme for bases.
	viridis_bases <- make_col_scheme(
		chars = c("A", "C", "G", "T"),
		groups = c("A", "C", "G", "T"),
		cols = c("#431352", "#34698c", "#44b57b", "#fde540")
	)

	## Make sequence logo.
	p <- consensus_matrix %>%
		pmap(function(name, value) {
			ggseqlogo(value, ncol = 1) +
				theme(text = element_text(size = font_size))
		})

	p <- plot_grid(
		plotlist = p,
		labels = pull(consensus_matrix, name),
		ncol = ncol,
		label_size = 5
	)

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

	## Start preparing data for plotting.
	seq_data <- tss_sequences %>%
		map(~ as.data.frame(.) %>% as_tibble(.name_repair = "unique", rownames = "name")) %>%
		bind_rows(.id = "sample") %>%
		rename(sequence = x)

	seq_data <- seq_data %>% 
		separate(
			name,
			into = c("chr", "start", "end", "strand", "score", "ntile"),
			sep = "_",
			remove = FALSE
		)

	## Get sequence length.
	sequence_length <- seq_data %>%
		pull(sequence) %>%
		nchar %>%
		unique / 2

	## Make columns for individual bases.
	base_columns <- seq_data %>%
		pull(sequence) %>%
		str_split("", simplify = TRUE)

	colnames(base_columns) <- as.character(c(-sequence_length:-1, 1:sequence_length))
	base_columns <- as_tibble(base_columns, .name_repair = "unique")
	seq_data <- bind_cols(seq_data, base_columns)

	## Get order of TSSs for plotting.
	tss_order <- seq_data %>%
		mutate(name = factor(name)) %>%
		group_by(chr, start, end, strand) %>%
		mutate(avg_score = mean(as.numeric(score))) %>%
		ungroup %>%
		mutate(name = fct_reorder(name, avg_score)) %>%
		pull(name) %>%
		levels

	## Format data for plotting.
	plot_data <- seq_data %>%
		gather(
			key = "position", value = "base",
			-sample, -name, -chr, -start, -end,
			-strand, -score, -ntile, -sequence
		)

	plot_data <- mutate(
		plot_data,
		position = as.numeric(position),
		name = factor(name, levels = tss_order),
		base = factor(base, levels = c("A", "C", "G", "T"))
	)

	n_samples <- plot_data %>%
		pull(sample) %>%
		unique %>%
		length

	## Plot sequence colormap
	p <- ggplot(plot_data, aes(x = position, y = name)) +
		geom_tile(aes(fill = base, color = base), ...) +
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
		)

	p <- p + facet_wrap(fct_rev(factor(ntile)) ~ sample, scales = "free", ncol = ncol)

	return(p)
}
