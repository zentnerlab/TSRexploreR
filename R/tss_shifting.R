#' TSS Shifting
#'
#' @description 
#' Analyze TSS shifts between samples within a consensus TSR set.
#'
#' @inheritParams common_params
#' @param sample_1 First sample to compare.
#'   Vector with sample name for TSS and TSR,
#'   with names 'TSS' and 'TSR'.
#' @param sample_2 Second sample to compare.
#'   Vector with sample name for TSS and TSR,
#'   with names 'TSS' and 'TSR'.
#' @param comparison_name Name assigned to the results in the TSRexploreR object.
#' @param tss_threshold Minimum number of raw counts required at a TSS for it to
#'   be considered in the shifting analysis.
#' @param max_distance TSRs less than this distance apart will be merged.
#' @param min_threshold Minimum number of raw counts required in each TSR for both samples.
#' @param n_resamples Number of resamplings for permutation test.
#'
#' @details
#' This function assesses the difference between TSS distributions from two distinct samples
#' in a set of consensus TSRs by calculating the earth mover's distance (EMD) between
#' them. For this approach, we imagine that the two TSS distributions in questions are piles 
#' of dirt, and ask how much dirt from one pile we would need to move, how far, and in which 
#' direction, to mimic the distribution of the other sample. The resulting score is between 
#' -1 and 1, with larger magnitudes indicating larger shifts and the sign indicating direction 
#' (negative values indicate upstream shifts and positive values indicate downstream shifts). 
#' The function also calculates a p-value for the null hypothesis that there is no difference 
#' (EMD = 0) based on a permutation test.
#'
#' 'sample_1' and 'sample_2' should be the names of the two samples to compare.
#' For the directions to make sense for the results 'sample_1' should be the control.
#'   and 'sample_2' the treatment sample.
#' The results will be stored back in the TSRexploreR object with the name given by
#'   'comparison_name'.
#' 'tss_threshold' applies a global threshold to remove TSSs below a certain score,
#'   and 'min_threshold' is the minimal score that both TSRs must have to be considered.
#' 'max_distance' is the maximum distance between two two TSRs to be considered for shifting.
#'
#' @return TSRexploreR object with shifting scores added.
#'
#' @examples
#' data(TSSs)
#'
#' tsre <- TSSs[c(1, 4)] %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3)
#' tss_shift(
#'   tsre,
#'   sample_1=c(TSS="S288C_WT_1", TSR="S288C_WT_1"),
#'   sample_2=c(TSS="S288C_D_1", TSR="S288C_D_1"),
#'   comparison_name="Untreated_vs_Diamide",
#'   max_distance = 100, min_threshold = 10, n_resamples = 1000L
#' )
#' 
#' @export

tss_shift <- function(
  experiment,
  sample_1,
  sample_2,
  comparison_name,
  tss_threshold=NULL,
  max_distance=100,
  min_threshold=10,
  n_resamples=1000L,
  fdr_cutoff=0.05
){

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    (is.character(sample_1) && length(sample_1) == 2) &&
    has_name(sample_1, c("TSS", "TSR"))
  )
  assert_that(
    (is.character(sample_2) && length(sample_2) == 2) &&
    has_name(sample_2, c("TSS", "TSR"))
  )
  assert_that(
    is.null(tss_threshold) ||
    (is.numeric(tss_threshold) && tss_threshold >= 0)
  )
  assert_that(is.count(max_distance))
  assert_that(is.count(min_threshold) && min_threshold > 5)
  assert_that(is.integer(n_resamples) && n_resamples >= 100L)
  assert_that(
    is.numeric(fdr_cutoff) &&
    (fdr_cutoff <= 1 & fdr_cutoff > 0)
  )
  assert_that(is.string(comparison_name))

  ## Retrieve TSSs and TSRs.
  TSSs <- extract_counts(experiment, "tss", c(sample_1["TSS"], sample_2["TSS"]))
  TSRs <- extract_counts(experiment, "tsr", c(sample_1["TSR"], sample_2["TSR"]))

  ## Merge TSRs into consensus ranges.
  consensus_TSRs <- TSRs %>%
    map(as_granges) %>%
    bind_ranges %>%
    GenomicRanges::reduce(min.gapwidth=max_distance, ignore.strand=FALSE) %>%
    as.data.table(key=c("seqnames", "strand", "start", "end"))
  consensus_TSRs[, FHASH := str_c(seqnames, start, end, strand, sep=":")]

  ## Remove consensus ranges width a width of 1.
  consensus_TSRs <- consensus_TSRs[start != end]

  ## Prepare TSSs for overlap with consensus ranges.
  TSSs <- rbindlist(TSSs, idcol="sample")
  setkey(TSSs, seqnames, strand, start, end)

  ## Filter out TSSs below threshold if requested.
  if (!is.null(tss_threshold)) {
    TSSs <- TSSs[score >= tss_threshold]
  }

  # Associate TSSs with consensus ranges.
  overlap <- foverlaps(TSSs, consensus_TSRs)
  overlap <- overlap[!is.na(start)]

  ## Filter out TSRs without TSSs in both samples.
  overlap[, count := uniqueN(sample), by=FHASH]
  overlap <- overlap[count == 2, .(sample, seqnames, start=i.start, end=i.end, strand, score, FHASH)]

  ## Get relative distances of each TSS in a TSR.
  overlap[, distance := ifelse(strand == "+", start - min(start), max(start) - start), by=FHASH]

  ## Prepare table for shift score calculation.
  overlap <- overlap %>%
    as_granges %>%
    sort %>%
    as.data.table
  overlap <- overlap[, .(sample, seqnames, start, end, strand, score, FHASH, distance)]

  ## Calculate the shift scores.
  shifts <- ShiftScores(
    overlap,
    baseline_level=sample_1["TSS"],
    nresamp=n_resamples,
    nthresh=min_threshold
  )

  ## p-value correction for multiple comparisons.
  setDT(shifts)
  shifts[, FDR := p.adjust(pval, "fdr")]
  shifts <- shifts[order(FDR)]

  ## Switch signs on shifting score so upstream is negative and
  ## downstream is positive.
  shifts[, shift_score := shift_score * -1]

  ## Filter out non-significant results.
  shifts <- shifts[FDR < fdr_cutoff]

  ## Add results to TSRexploreR object.
  shifts[, c("start", "end") := list(as.numeric(start), as.numeric(end))]
  experiment@shifting$results[[comparison_name]] <- shifts

  return(experiment)
}

#' Shift Score
#'
#' Calculate shift scores and associated permutation test p-values.
#'
#' @importFrom Rcpp sourceCpp
#'
#' @param tss_table Table of TSSs prepared for shifting score calculation.
#' @param baseline_level The sample being used as the baseline for calculation.
#' @param calc_pvalue Whether p-values should be returned for comparisons.
#' @param nresamp Number of resamplings for permutation test.
#' @param nthresh Minimum number of raw counts required in each TSR for both samples.
#' @param check_sort Check that the input is sorted properly.
#'
#' @rdname ShiftScores-function
#' @export

ShiftScores <- function(
  tss_table,
  baseline_level,
  calc_pvalue=TRUE, 
  nresamp=100L,
  nthresh=2,
  check_sort=TRUE
){

  dat <- tss_table %>%
    as_tibble %>%
    dplyr::select(fhash=FHASH, sample_indicator=sample, distances=distance, scores=score)
  
  # Assumes fhash is consecutive, no regrouping necessary.
  # Assumes there are only two samples.
  dat <- dplyr::mutate(
    dat, sample_indicator=as.integer(sample_indicator==baseline_level)
  )

  if(check_sort) dat <- dplyr::arrange(dat, fhash, sample_indicator, distances)

  ## Filter out TSRs where one of the samples has too low a score.
  out_frame <- dat %>%
    dplyr::group_by(fhash, sample_indicator) %>%
    dplyr::summarise(sum_score=sum(scores)) %>% 
    dplyr::mutate(toosmall= sum_score < nthresh)

  if(sum(out_frame$toosmall) > 0){
    warning("Some sequences have fewer than nthresh scores for at least one sample. 
            These are ignored and returned as NA.")
  }

  dat <- dat %>%
    dplyr::left_join(out_frame) %>%
    dplyr::group_by(fhash) %>%
    dplyr::filter(!any(toosmall)) %>%
    dplyr::ungroup()

  ## Filter out TSRs where only one sample is present.
  dat <- dat %>%
    dplyr::group_by(fhash) %>%
    dplyr::filter(dplyr::n_distinct(sample_indicator) == 2) %>%
    dplyr::ungroup()

  ## Calculate the shift score.
  out <- with(
    dat, ## returns a 2 by n_distinct matrix
    allTheShiftScores(
      fhash, distances, scores, sample_indicator, 
      as.integer(calc_pvalue), as.integer(nresamp), dplyr::n_distinct(fhash)
    )
  )

  ## Add the coordinates back to the shift score.
  
  ## Altered to avoid warning message on naming.
  out <- t(out)
  colnames(out) <- c("shift_score", "pval")
  outdf <- as_tibble(out)
  # outdf <- out %>%
  #   t %>%
  #   as_tibble(.name_repair="unique") %>%
  #   dplyr::rename(shift_score=1, pval=2)

  outdf <- dat %>%
     dplyr::distinct(fhash) %>%
     dplyr::bind_cols(outdf) %>%
     tidyr::separate(fhash, into=c("seqnames", "start", "end", "strand"), sep=":")

  return(outdf)
}
