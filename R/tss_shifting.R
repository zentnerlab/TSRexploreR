
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

	## Calculate the shift scores.
	shifts <- ShiftScores(
		overlap$FHASH, overlap$sample, overlap$score, nresamp = 1000L,
		baseline_level = "Untreated", nthresh = 10
	)
}

#' Shifting Score
#'
#' Calculating shifting score and p-value.
#'
#' @import RcppArmadillo
#'
#' @param fhash fhash of set
#' @param sample_indicator The column with the two sample names
#' @param distances bin positions
#' @param scores bin scores
#' @param calc_pvalue Should p-value be returned for comparison
#' @param nresamp number of resamples
#' @param baseline_level control smaple
#' @param nthresh both samples must have at least this number of reads
#' @param check_sort Check that the input is sorted properly
#'
#' @rdname ShiftScores-function
#' @export

ShiftScores <- function(
	fhash, sample_indicator, distances, scores, calc_pvalue = TRUE, 
	nresamp = 100L, baseline_level = sample_indicator[1], nthresh = 2,
	check_sort = TRUE
){

  dat = data.frame(fhash, sample_indicator, distances, scores)
  if(check_sort) dat = dplyr::arrange(dat, fhash, sample_indicator, distances)
  
  # Assumes fhash is consequtive, no regrouping necessary
  # Assumes there are only two samples
  dat = dplyr::mutate(
    dat, sample_indicator = as.integer(sample_indicator==baseline_level)
    )
  out_frame = dat %>% dplyr::group_by(fhash,sample_indicator) %>%
    dplyr::summarise(n = dplyr::n()) %>% 
    dplyr::summarise(smallest = min(n)) %>%
    dplyr::mutate(toosmall = smallest < nthresh)
  if(sum(out_frame$toosmall) > 0){
    warning("Some sequences have fewer than nthresh scores for at least one sample. 
            These are ignored and returned as NA.")
  }
  dat = dplyr::left_join(dat, out_frame) %>% filter(!toosmall)
  out = with(dat, ## returns a 2 by n_distinct matrix
    allTheShiftScores(fhash, distances, scores, sample_indicator, 
                      as.integer(calc_pvalue), nresamp, dplyr::n_distinct(fhash))
  )
  outdf = data.frame(shift_score = out[1,])
  if(calc_pvalue) outdf$pval = out[2,]
  outdf = dplyr::bind_cols(
    outdf, out_frame %>% dplyr::filter(!toosmall) %>% dplyr::select(fhash)
    )
  outdf = out_frame %>% dplyr::select(fhash) %>% dplyr::left_join(outdf) 
  outdf
}

### Some tests
# big_dat = read_delim("../../shifting_score.tsv", delim = "\t")
#real_small_dat = big_dat %>% filter(FHASH == "b0cdc3a9f9e98f26b0325802c184a975") %>%
#  arrange(sample, distance)
#
#small_but_ugly = big_dat[1:67,]
#med_dat = big_dat[1:2028,]
#
#with(real_small_dat,
#  allTheShiftScores(FHASH, distance, score, as.integer(sample=="Untreated"), 1, 10, 1))
#
#with(real_small_dat, ShiftScores(FHASH, sample, distance, score))
#
#with(small_but_ugly, ShiftScores(FHASH, sample, distance, score))
#
#system.time(with(med_dat, ShiftScores(FHASH, sample, distance, score, nresamp = 1000)))
#
#biggish_dat = big_dat[1:19997,]
#system.time(out <- with(biggish_dat, ShiftScores(FHASH, sample, distance, score, nresamp = 1000)))
#}
