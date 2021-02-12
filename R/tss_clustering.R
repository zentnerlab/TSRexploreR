#' TSS Clustering
#'
#' @description
#' Basic distance and threshold-based clustering of TSSs.
#'
#' @inheritParams common_params
#' @param max_distance Maximum allowable distance between TSSs for clustering.
#' @param max_width Maximum allowable TSR width.
#' @param n_samples Keep TSS if 'threshold' number of counts is present in n_samples
#'   number of samples.
#' @param singlet_threshold TSRs of width 1 must have a score greater than
#'   or equal to this threshold to be kept.
#'
#' @details
#' Genes rarely have a single TSS, but rather a cluster of TSSs.
#' This function clusters TSSs into Transcription Start Regions (TSRs).
#' TSSs are clustered if their score is greater than or equal to 'threshold'
#'   in at least 'n_samples' number of samples,
#'   and are less than or equal to 'max_distance' from each other.
#' The clustered TSSs cannot encompass more than 'max_width' bases.
#'   A global singlet threshold can be applied using 'singlet_threshold'.
#'
#' @return TSRexploreR object with TSRs added to GRanges and data.table counts.
#'
#' @examples
#' data(TSSs_reduced)
#'
#' exp <- TSSs_reduced %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss")
#' exp <- tss_clustering(exp, threshold=3)
#'
#' @export

tss_clustering <- function(
  experiment,
  samples="all",
  threshold=NULL,
  n_samples=NULL,
  max_distance=25,
  max_width=NULL,
  singlet_threshold=NULL
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.count(max_distance))
  assert_that(is.null(max_width) || is.count(max_width))
  assert_that(is.null(n_samples) || is.count(n_samples))
  assert_that(
    is.null(singlet_threshold) ||
    (is.numeric(singlet_threshold) && singlet_threshold > 0)
  )

  ## Get selected samples.
  select_samples <- extract_counts(experiment, "tss", samples, FALSE)

  ## Filter samples if requested.
  if (!is.null(threshold)) {
    if (!is.null(n_samples)) {
      select_samples <- rbindlist(select_samples, idcol="samples")
      select_samples[, .n_samples := .N, by=FHASH]
      select_samples <- select_samples[score >= threshold & .n_samples >= n_samples]
      select_samples[, .n_samples := NULL]
      select_samples <- split(select_samples, by="samples", keep.by=FALSE)
    } else {
      select_samples <- map(select_samples, ~.x[score > threshold])
    }
  }

  select_samples <- map(select_samples, function(x) {
    keep_cols <- c("seqnames", "start", "end", "strand", "score")
    if (any(colnames(x) == "normalized_score")) {
      keep_cols <- c(keep_cols, "normalized_score")
    }
    x <- x[, ..keep_cols]
    x <- as_granges(x)
    return(x)
  })
  
  ## Cluster TSSs into TSRs.
  clustered_TSSs <- map(
    select_samples, .aggr_scores,
    max_distance, max_width, singlet_threshold
  )

  ## Add TSRs back to TSRexploreR object.
  clustered_TSSs <- map(clustered_TSSs, function(x) {
    x <- as.data.table(x)
    x[,
      FHASH := str_c(seqnames, start, end, strand, sep=":"),
      by=seq_len(nrow(x))
    ]
    return(x)
  })

  TSR_granges <- map(clustered_TSSs, function(x) {
    keep_cols <- c("seqnames", "start", "end", "strand", "score")
    if (any(colnames(x) == "normalized_score")) {
      keep_cols <- c(keep_cols, "normalized_score")
    }
    TSR_granges <- x[, ..keep_cols]
    TSR_granges <- as_granges(TSR_granges)
    return(TSR_granges)
  })

  experiment@experiment$TSRs <- TSR_granges
  experiment <- set_count_slot(
    experiment, clustered_TSSs,
    "counts", "tsr", "raw"
  )
  return(experiment)

}

## Aggregate Scores
##
## @param granges GRanges.
## @param maxdist Maximum distance to cluster.
## @param sthresh Singlet threshold.
## @param maxwidth Maximum allowable TSR width.

.aggr_scores <- function(granges, maxdist, maxwidth, sthresh) {

  ## Check inputs.
  assert_that(is(granges, "GRanges"))
  assert_that(is.count(maxdist) | maxdist == 0)

  ## Get TSSs within max distance.
  reduced <- granges %>%
    stretch(maxdist * 2) %>%
    reduce_ranges_directed

  granges <- as.data.table(granges, key=c("seqnames", "strand", "start", "end"))
  reduced <- as.data.table(reduced, key=c("seqnames", "strand", "start", "end"))

  overlaps <- foverlaps(reduced, granges)

  ## Calculate aggregate score and number of unique TSSs per TSR.
  if ("normalized_score" %in% colnames(overlaps)) {
    overlaps <- overlaps[,
      .(
        start=min(start), end=max(end), score=sum(score), n_unique=.N,
        normalized_score=sum(normalized_score)
      ),
      by=.(seqnames, strand, i.start, i.end)
    ]
  } else {
    overlaps <- overlaps[,
      .(start=min(start), end=max(end), score=sum(score), n_unique=.N),
      by=.(seqnames, strand, i.start, i.end)
    ]
  }
  overlaps[, c("i.start", "i.end") := NULL]
  overlaps <- overlaps %>%
    as_granges %>%
    sort %>%
    as.data.table

  ## If max_width is set, remove TSRs that are too wide.
  if (!is.null(maxwidth)) {
    overlaps <- overlaps[width <= maxwidth]
  }

  ## If singlet threshold is set,
  ## remove singlets that fail threshold check.
  if (!is.null(sthresh)) {
    overlaps <- overlaps[!(width == 1 & score < sthresh)]
  }

  return(overlaps)

}
