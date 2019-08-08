
#' Retrieve Sequences Near TSSs
#'
#' Retrieve sequences surrounding TSSs for further plotting
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom Rsamtools FaFile
#' @importFrom Biostrings getSeq DNAStringSet
#' @importFrom GenomicRanges GRanges resize score width seqnames start strand
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param sample Name of sample to analyze
#' @param genome_assembly Genome assembly in fasta format
#' @param distance Bases to add on each side of TSS
#'
#' @return Sequences surrounding TSSs
#'
#' @export
#' @rdname tss_sequences-function

tss_sequences <- function(experiment, sample, genome_assembly, threshold = 1, distance = 10) {
	## Open genome assembly
	assembly <- FaFile(genome_assembly)

	## Get TSS sequences.
	tss_sequences <- experiment@experiment$TSSs[[sample]] %>%
		.[score(.) >= threshold]

	tss_sequence_names <- paste(
		as.character(seqnames(tss_sequences)), 
		start(tss_sequences),
		strand(tss_sequences),
		score(tss_sequences),
		sep = "_"
	)
		
	tss_sequences <- tss_sequences %>%
		resize(width = distance, fix = "start") %>%
		resize(width = width(.) + distance, fix = "end") %>%
		getSeq(assembly, .)

	names(tss_sequences) <- tss_sequence_names

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
#'
#' @param tss_sequences Sequences surrounding TSS generated with tss_sequences
#'
#' @return ggplot2 object with sequence logo
#'
#' @export
#' @rdname plot_sequence_logo-function

plot_sequence_logo <- function(tss_sequences) {
	## Calculate consensus matrix.
	consensus_matrix <- tss_sequences %>%
		consensusMatrix(.,as.prob = TRUE) %>%
		.[rownames(.) %in% c("A", "C", "G", "T"), ]

	## Create viridis color scheme for bases.
	viridis_bases <- make_col_scheme(
		chars = c("A", "C", "G", "T"),
		groups = c("A", "C", "G", "T"),
		cols = c("#431352", "#34698c", "#44b57b", "#fde540")
	)

	## Make sequence logo.
	p <- ggplot() +
		geom_logo(consensus_matrix, col_scheme = viridis_bases) +
		theme_logo()

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
	p <- ggplot(tss_sequences, aes(x=Position, y=sequence)) +
		geom_tile(aes(fill=base)) +
		scale_fill_viridis_d() +
		theme_minimal() +
		theme(
			axis.title.y=element_blank(),
			axis.text.y=element_blank(),
			legend.title=element_blank(),
			axis.title.x=element_text(size=16, margin=margin(t=15)),
			panel.grid=element_blank()
		)

	return(p)
}
