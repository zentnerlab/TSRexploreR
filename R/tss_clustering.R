
#' TSS Clustering
#'
#' @description
#' Basic distance and threshold-based clustering of TSSs.
#'
#' @param experiment tsrexplorer object
#' @param threshold Consider only TSSs with at least this number of raw counts
#' @param samples Samples for which TSSs should be clustered
#' @param max_distance Maximum distance between TSSs that can be clustered
#'
#' @details
#' This function clusters TSSs into Transcription Start Regions (TSRs).
#'
#' TSSs are clustered if their score is greater than or equal to
#'   'threshold' and are less than or equal to 'max_distance'
#'   from each other.
#'
#' @return tsrexplorer object with TSRs
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#'
#' @rdname tss_clustering-function
#' @export

tss_clustering <- function(
  experiment,
  samples = "all",
  threshold = 1,
  max_distance = 25
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(is.count(threshold) && threshold > 0)
  assert_that(is.count(max_distance) && max_distance >= 5)

  ## Get selected samples.
  select_samples <- extract_counts(experiment, "tss", samples)
  
  select_samples <- select_samples %>%
    map(function(x) {
      x <- x[
        score >= threshold,
        .(seqnames, start, end, strand, score)
      ]
      return(x)
    })

  ## Convert samples to GRanges.
  select_samples <- map(select_samples, as_granges)

  ## Call TSRs.
  clustered_TSSs <- select_samples %>%
    map(function(x) {

      # Cluster TSSs within 'max_distance'.
      clustered <- GenomicRanges::reduce(
        x, with.revmap = TRUE,
        min.gapwidth = max_distance + 1
      )

      # Get aggregate sum of scores.
      cluster_info <- aggregate(
        x, mcols(clustered)$revmap,
        score = sum(score),
        n_unique = lengths(score)
      )

      clustered$score <- cluster_info$score
      clustered$n_unique <- cluster_info$n_unique
      clustered$revmap <- NULL

      return(clustered)
    })

  ## Add TSRs back to tsrexplorer object.
  clustered_TSSs <- map(clustered_TSSs, function(x) {
    x <- as.data.table(x)
    x[, FID := seq_len(nrow(x))]
    x[,
      FHASH := str_c(seqnames, start, end, strand, sep = ":"),
      by = seq_len(nrow(x))
    ]
    return(x)
  })

  TSR_granges <- map(clustered_TSSs, function(x) {
    TSR_granges <- x[, .(seqnames, start, end, strand, score)]
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
