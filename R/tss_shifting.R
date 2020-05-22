
#' TSS Shifting
#'
#' Calculate TSS shifting statistics
#'
#' @param experiment tsrexplorer object
#' @param compare_samples Vector of names of the two TSR samples to compare
#' @param min_distance TSRs less than this distance apart will be merged
#' @param min_threshold Minimum number of raw counts required in each TSR for both TSR samples
#' @param n_resamples Number of resamplings for permutation test
#'
#' @rdname tss_shift-function
#' @export

tss_shift <- function(
	experiment,
	compare_samples,
	min_distance = 100,
	min_threshold = 10,
	n_resamples = 1000L
){

	## Input checks.
	if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")

	if (!is(compare_samples, "character") || length(compare_samples) != 2) {
		stop("compare_samples must be a character vector with the names of the two TSR samples to compare")
	}

	if (
		!is(min_distance, "numeric") || length(min_distance) > 1 ||
		min_distance %% 1 != 0 || min_distance < 0
	) {
		stop("min_distance must be 0 or a positive integer")
	}

	if (
                !is(min_threshold, "numeric") || length(min_threshold) > 1 ||
                min_threshold %% 1 != 0 || min_threshold < 1
	) {
		stop("min_threshold must be a positive integer")
	}

	if (!is(n_resamples, "integer") || length(n_resamples) > 1 || n_resamples < 100) {
		stop("n_resamples must be a positive integer greater than or equal to 100")
	}

	
	## Get samples to be compared.
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

	## Associate consensus TSRs with TSSs.
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

	## Filter out TSRs without TSSs in both samples.
	overlap[, count := uniqueN(sample), by = FHASH]
	overlap <- overlap[count == 2, .(sample, seqnames, start, end, strand, score, FHASH)]

	## Get relative distances of each TSS in a TSR.
	overlap[, distance := ifelse(strand == "+", start - min(start), max(start) - start), by = FHASH]

	## Calculate the shift scores.
	shifts <- ShiftScores(
		fhash = overlap$FHASH, sample_indicator = overlap$sample, 
		distances = overlap$distance, scores = overlap$score,
		nresamp = n_resamples, baseline_level = compare_samples[1],
		nthresh = min_threshold
	)

	shifts <- as.data.table(shifts)
	setnames(shifts, "fhash", "FHASH")
	shifts[, FDR := p.adjust(pval, "fdr")]

	## Merge data back into the tsrexplorer object.
	shift_results <- as.data.table(consensus_tsrs)
	shift_results <- merge(shift_results, shifts, by = "FHASH")

	## Return table of results.
	shift_results[, FHASH := NULL]
	shift_results <- shift_results[order(FDR)]

	return(shift_results)
}

#' Shifting Score
#'
#' Calculate shifting scores and associated permutation test p-values.
#'
#' @importFrom Rcpp sourceCpp
#'
#' @param fhash fhash of set
#' @param sample_indicator The column with the names of the samples to be compared (where is the column from qq)
#' @param distances bin positions
#' @param scores bin scores
#' @param calc_pvalue Whether p-values should be returned for comparisons
#' @param nresamp Number of resamplings for the permutation test
#' @param baseline_level Control smaple
#' @param nthresh Both samples must have at least this number of reads in each TSR
#' @param check_sort Check that the input is sorted properly (by fhash? qq)
#'
#' @rdname ShiftScores-function
#' @export

ShiftScores <- function(
  fhash, sample_indicator, distances, scores, calc_pvalue = TRUE, 
  nresamp = 100L, baseline_level = sample_indicator[1], nthresh = 2,
  check_sort = TRUE){
  
  dat = data.frame(fhash, sample_indicator, distances, scores)
  if(check_sort) dat = dplyr::arrange(dat, fhash, sample_indicator, distances)
  
  # Assumes fhash is consecutive, no regrouping necessary.
  # Assumes there are only two samples.
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
                               as.integer(calc_pvalue), as.integer(nresamp), dplyr::n_distinct(fhash))
  )
  outdf = data.frame(shift_score = out[1,])
  if(calc_pvalue) outdf$pval = out[2,]
  outdf = dplyr::bind_cols(
    outdf, out_frame %>% dplyr::filter(!toosmall) %>% dplyr::select(fhash)
  )
  outdf = out_frame %>% dplyr::select(fhash) %>% dplyr::left_join(outdf) 
  return(outdf)
}
