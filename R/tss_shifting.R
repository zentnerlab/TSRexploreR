
#' TSS Shifting
#'
#' Calculate TSS shifting statistics
#'
#' @param experiment tsrexplorer object
#' @param compare_samples Vector of two samples to compare
#'
#' @rdname tss_shift-function
#' @export

tss_shift <- function(experiment, compare_samples, min_distance = 100) {
	
	## Grab samples to be compared.
	select_samples <- extract_counts(experiment, "tss", compare_samples)

	## Get consensus TSRs.
	consensus_tsrs <- rbindlist(select_samples)[, .(tsr_coords)]
	consensus_tsrs <- unique(consensus_tsrs)
	consensus_tsrs[,
		c("seqnames", "start", "end", "strand") :=
		tstrsplit(tsr_coords, split = ":")
	]

	consensus_tsrs <- makeGRangesFromDataFrame(consensus_tsrs)
	consensus_tsrs <- GenomicRanges::reduce(
		consensus_tsrs, ignore.strand = FALSE, min.gapwidth = min_distance
	)
	
	consensus_tsrs <- as.data.table(consensus_tsrs)
	consensus_tsrs[,
		FHASH := digest(str_c(seqnames, start, end, strand, collapse = ":")),
		by = seq_len(nrow(consensus_tsrs))
	]
	consensus_tsrs <- makeGRangesFromDataFrame(consensus_tsrs, keep.extra.columns = TRUE)

	## Assign consensus TSR to TSSs.
	tss_data <- rbindlist(select_samples, idcol = "sample")[,
		.(sample, seqnames, start, end, strand, score)
	]
	tss_data <- makeGRangesFromDataFrame(tss_data, keep.extra.columns = TRUE)

	overlap <- findOverlapPairs(query = consensus_tsrs, subject = tss_data)
	overlap <- as.data.table(overlap)[,
		.(second.sample, second.X.seqnames, second.X.start, second.X.end,
		second.X.strand, second.X.score, first.X.FHASH)
	]

	setnames(
		overlap,
		old = c(
			"second.sample", "second.X.seqnames", "second.X.start",
			"second.X.end", "second.X.strand", "second.X.score",
			"first.X.FHASH"
		),
		new = c("sample", "seqnames", "start", "end", "strand", "score", "FHASH")
	)

	## Filter out TSRs that don't have TSSs in both samples.
	overlap[, count := uniqueN(sample), by = FHASH]
	overlap <- overlap[count == 2, .(sample, seqnames, start, end, strand, score, FHASH)]

	## Get relative distances of each TSS in a TSR.
	overlap[, distance := ifelse(strand == "+", start - min(start), max(start) - start), by = FHASH]
}
