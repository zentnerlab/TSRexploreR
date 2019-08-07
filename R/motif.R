
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
	tss_sequences <- experiment@experiment[[sample]] %>%
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
