
#' Retrieve Sequences Near TSSs
#'
#' Retrieve sequences surrounding TSSs for further plotting
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom Rsamtools FaFile
#' @importFrom purrr map map2
#' @importFrom dplyr mutate filter ntile pull
#' @importFrom Biostrings getSeq DNAStringSet
#' @importFrom GenomicRanges GRanges resize score width seqnames start strand makeGRangesFromDataFrame
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param sample Either "all" or names of samples to analyze
#' @param genome_assembly Genome assembly in fasta format
#' @param distance Bases to add on each side of TSS
#' @param quantiles Break data into quantiles
#'
#' @return Sequences surrounding TSSs
#'
#' @export
#' @rdname tss_sequences-function

tss_sequences <- function(experiment, samples = "all", genome_assembly, threshold = 1, distance = 10, quantiles = 1) {

	## Open genome assembly.
	genome_assembly <- FaFile(genome_assembly)

	## Pull selected samples.
	if (samples == "all") samples <- names(experiment@experiment$TSSs)

	## Start preparing data for retrieving sequences.
	tss_sequences <- experiment@experiment$TSSs[samples] %>%
		map(~ as_tibble(., .name_repair = "unique")) %>%
		bind_rows(.id = "sample") %>%
		filter(score >= threshold)

	if (quantiles > 1) tss_sequences <- mutate(tss_sequences, ntile = ntile(score, quantiles))

	## Store sequence names.
	tss_sequence_names <- tss_sequences %>%
		split(.$sample) %>%
		map(~ mutate(., tss_name = paste(seqnames, start, end, width, strand, score, sep = "_")))

	if (quantiles > 1) {
		tss_sequence_names <- map(
			tss_sequence_names,
			~ mutate(., tss_name = paste(tss_name, ntile, sep = "_"))
		)
	}

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

	## Return result.
	if (quantiles > 1) {
		tss_sequences <- list(
			quantile_plot = TRUE,
			tss_sequences = tss_sequences
		)
	} else {
		tss_sequences <- list(
			quantile_plot = FALSE,
			tss_sequences = tss_sequences
		)
	}

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
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @param tss_sequences Sequences surrounding TSS generated with tss_sequences
#' @param ncol Number of columns to plot if quantiles is not set
#'
#' @return ggplot2 object with sequence logo
#'
#' @export
#' @rdname plot_sequence_logo-function

plot_sequence_logo <- function(tss_sequences, ncol = 1) {

	## Grab sequences from input.
	tss_seqs <- tss_sequences$tss_sequences

	## Calculate consensus matrix.
	consensus_matrix <- tss_seqs %>%
		map(
			~ consensusMatrix(., as.prob = TRUE) %>%
			.[rownames(.) %in% c("A", "C", "G", "T"), ]
		)

	## Create viridis color scheme for bases.
	viridis_bases <- make_col_scheme(
		chars = c("A", "C", "G", "T"),
		groups = c("A", "C", "G", "T"),
		cols = c("#431352", "#34698c", "#44b57b", "#fde540")
	)

	## Make sequence logo.
	p <- ggplot() +
		geom_logo(consensus_matrix) +
		theme_logo() +
		facet_wrap(~ seq_group, ncol = ncol, scales = "free_x")

	return(p)
}

#' Plot Sequence Colormap
#'
#' Make a sequence colormap for the sequences around TSSs.
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom Biostrings DNAStringSet
#' @importFrom BiocGenerics width
#'
#' @param tss_sequences Sequences surrounding TSS generated with tss_sequences
#'
#' @return ggplot2 object of sequence colormap
#'
#' @export
#' @rdname plot_sequence_colormap-function

plot_sequence_colormap <- function(tss_sequences) {
	## Get sequence length
	sequence_length <- tss_sequences %>%
		width(.) %>% unique / 2

	## Format data for plotting
	tss_sequences <- tss_sequences %>%
		as.character %>%
		str_split(pattern = "", simplify = TRUE) %>%
		as_tibble(.name_repair="unique") %>%
		setNames(c(-sequence_length:-1, 1:sequence_length)) %>%
		rowid_to_column(var = "sequence") %>%
		gather(key = "Position", value = "base", -sequence) %>%
		mutate(Position = factor(Position, levels=c(-sequence_length:-1, 1:sequence_length)))

	## Plot sequence colormap
	p <- ggplot(tss_sequences, aes(x = Position, y = sequence)) +
		geom_tile(aes(fill = base)) +
		scale_fill_viridis_d() +
		theme_minimal() +
		theme(
			axis.title.y=element_blank(),
			axis.text.y=element_blank(),
			legend.title=element_blank(),
			axis.title.x=element_text(margin = margin(t = 15)),
			panel.grid=element_blank()
		)

	return(p)
}
