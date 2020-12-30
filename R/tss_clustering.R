#' TSS Clustering
#'
#' @description
#' Basic distance and threshold-based clustering of TSSs.
#'
#' @inheritParams common_params
#' @param max_distance Maximum distance between TSSs that can be clustered
#' @param max_width Maximum width of TSR allowed.
#'
#' @details
#' This function clusters TSSs into Transcription Start Regions (TSRs).
#'
#' TSSs are clustered if their score is greater than or equal to
#'   'threshold' and are less than or equal to 'max_distance'
#'   from each other.
#'
#' @return TSRexploreR object with TSRs
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#'
#' @rdname tss_clustering-function
#' @export

tss_clustering <- function(
  experiment,
  samples="all",
  threshold=NULL,
  max_distance=25,
  max_width=NULL
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.count(max_distance))
  assert_that(is.null(max_width) || is.count(max_width))

  ## Get selected samples.
  select_samples <- experiment %>%
    extract_counts("tss", samples, FALSE) %>%
    preliminary_filter(FALSE, threshold)

  select_samples <- map(select_samples, function(x) {
    keep_cols <- c("seqnames", "start", "end", "strand", "score")
    if (any(colnames(x) == "normalized_score")) {
      keep_cols <- c(keep_cols, "normalized_score")
    }
    x <- x[, ..keep_cols]
    x <- as_granges(x)
    return(x)
  })
  
  ## Call TSRs.
  clustered_TSSs <- map(select_samples, .aggr_scores, max_distance, max_width)

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

#' Aggregate Scores
#'
#' @param granges GRanges
#' @param maxdist Maximum distance to cluster

.aggr_scores <- function(granges, maxdist, maxwidth) {

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

  ## If max_width is set remove TSRs that are too wide.
  if (!is.null(maxwidth)) {
    overlaps <- overlaps[width <= maxwidth]
  }

  return(overlaps)

#  ## Cluster TSSs within 'max_distance'
#  clustered <- GenomicRanges::reduce(
#    granges, with.revmap=TRUE,
#    min.gapwidth=maxdist + 1
#  )
#
#  ## Get aggregate sum of scores.
#  if (any(colnames(mcols(granges)) == "normalized_score")) {
#    cluster_info <- aggregate(
#      granges, mcols(clustered)$revmap, 
#      score=sum(score),
#      normalized_score=sum(normalized_score),
#      n_unique=length(score)
#    )
#  } else {
#    cluster_info <- aggregate(
#      granges, mcols(clustered)$revmap,
#      score=sum(score),
#      n_unique=length(score)
#    )
#  }
#
#  clustered$score <- cluster_info$score
#  clustered$n_unique <- cluster_info$n_unique
#  if (any(colnames(mcols(granges)) == "normalized_score")) {
#    clustered$normalized_score <- cluster_info$normalized_score
#  }
#  clustered$revmap <- NULL
}
