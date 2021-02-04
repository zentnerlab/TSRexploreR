
#' Shifting Rank Plot
#'
#' @description
#' Barplot of shift ranks ordered by score.
#'
#' @inheritParams common_params
#' @param score_order Either descending or ascending
#'
#' @details
#' The 'tss_shifting' function uses the Earth Movers Distance (EMD)
#'   to find significant upstream or downstream shifts in TSS clusters.
#' EMD returns a score between -1 and 1 inclusive, with negative and positive scores
#'   correspending to upstream and downstream shifts respectively.
#' In rarer cases you can have a shifting scores of 0,
#'   indicating an equal shift of TSS cluster 'mass' upstream and downstream.
#'
#' This functions generates a barplot of shifting score,
#'   with the x-axis corresponding to TSS clusters ordered by tss shifting score,
#'   and the y-axis the shifting score.
#' TSS clusters can be in ascending or descending order as controlled by 'score_order'.
#'
#' This function 
#'
#' @return ggplot2 object of shifting rank plot.
#'
#' @seealso
#' \code{\link{tss_shift}} to calculate TSS cluster shifting.
#'
#' @examples
#' data(TSSs)
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "TSRexploreR")
#' samples <- data.frame(
#'   sample_name=c(sprintf("S288C_D_%s", seq_len(3)), sprintf("S288C_WT_%s", seq_len(3))),
#'   file_1=rep(NA, 6), file_2=rep(NA, 6),
#'   condition=c(rep("Diamide", 3), rep("Untreated", 3))
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=samples, genome_assembly=assembly) %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3) %>%
#'   merge_samples(data_type = "tss", merge_group="condition") %>%
#'   merge_samples(data_type = "tsr", merge_group="condition") %>%
#'   tss_shift(
#'     sample_1=c(TSS="S288C_WT_1", TSR="S288C_WT_1"),
#'     sample_2=c(TSS="S288C_D_1", TSR="S288C_D_1"),
#'     comparison_name="Untreated_vs_Diamide",
#'     max_distance = 100, min_threshold = 10, n_resamples = 1000L
#'   )
#' p <- plot_shift_rank(tsre)
#'
#' @export

plot_shift_rank <- function(
  experiment,
  samples="all",
  score_order="descending",
  ncol=3
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  score_order <- match.arg(
    str_to_lower(score_order),
    c("descending", "ascending")
  )
  assert_that(is.count(ncol))

  ## Retrieve samples.
  if (samples == "all") {
    select_samples <- experiment@shifting$results
  } else {
    select_samples <- experiment@shifting$results[samples]
  }

  select_samples <- rbindlist(select_samples, idcol="comparison")

  ## Order features by selected order.
  select_samples[, FID := as.character(seq_len(nrow(select_samples)))]

  if (score_order == "descending") {
    select_samples[, FID := fct_reorder(FID, shift_score, .desc=TRUE)]
  } else {
    select_samples[, FID := fct_reorder(FID, shift_score, .desc=FALSE)]
  }

  ## Set sample order if required.
  if (!all(samples == "all")) {
    select_samples[, comparison := factor(comparison, levels=samples)]
  }

  ## Generate the plot.
  p <- ggplot(select_samples, aes(x=.data$FID, y=.data$shift_score, fill=.data$shift_score)) +
    geom_col(width=1) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    ) +
    facet_wrap(~ comparison, ncol=ncol, scales="free")

  return(p)
}

#' Shift Count Plot
#'
#' @description
#' Generate a stacked barplot for the number of TSS clusters
#'   shifted upstream and/or downstream per sample.
#'
#' @inheritParams common_params
#' @param ... Arguments passed to geom_col
#'
#' @return ggplot2 of stacked barplot of shifted TSS clusters.
#'   If 'return_table' is TRUE, a data.frame of underlying counts
#'   is returned instead.
#'
#' @seealso
#' \code{\link{tss_shift}} For TSS cluster shifting calculation.
#'
#' @details
#' The 'tss_shifting' function uses the Earth Movers Distance (EMD)
#'   to find significant upstream or downstream shifts in TSS clusters.
#' EMD returns a score between -1 and 1 inclusive, with negative and positive scores
#'   correspending to upstream and downstream shifts respectively.
#' In rarer cases you can have a shifting scores of 0,
#'   indicating an equal shift of TSS cluster 'mass' upstream and downstream.
#'
#' This function generates a stacked barplot for the number of TSS clusters
#'   with an upstream and downstream shift.
#' If 'return_table' is TRUE, a data.frame is returned that provides the underlying
#'   counts for each sample.
#'
#' @examples
#' data(TSSs)
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "TSRexploreR")
#' samples <- data.frame(
#'   sample_name=c(sprintf("S288C_D_%s", seq_len(3)), sprintf("S288C_WT_%s", seq_len(3))),
#'   file_1=rep(NA, 6), file_2=rep(NA, 6),
#'   condition=c(rep("Diamide", 3), rep("Untreated", 3))
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=samples, genome_assembly=assembly) %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3) %>%
#'   merge_samples(data_type = "tss", merge_group="condition") %>%
#'   merge_samples(data_type = "tsr", merge_group="condition") %>%
#'   tss_shift(
#'     sample_1=c(TSS="S288C_WT_1", TSR="S288C_WT_1"),
#'     sample_2=c(TSS="S288C_D_1", TSR="S288C_D_1"),
#'     comparison_name="Untreated_vs_Diamide",
#'     max_distance = 100, min_threshold = 10, n_resamples = 1000L
#'   )
#' p <- plot_shift_count(tsre)
#'
#' @export

plot_shift_count <- function(
  experiment,
  samples="all",
  ncol=3,
  return_table=FALSE,
  ...
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(is.count(ncol))
  assert_that(is.flag(return_table))

  ## Get samples.
  if (all(samples == "all")) {
    select_samples <- experiment@shifting$results
  } else {
    select_samples <- experiment@shifting$results[samples]
  }

  select_samples <- rbindlist(select_samples, idcol="comparison")

  ## Annotate shifting status.
  select_samples[, shift_status := case_when(
    shift_score < 0 ~ "upstream",
    shift_score > 0 ~ "downstream",
    TRUE ~ "n.s."
  )]

  ## Get number of shifts per shifting status.
  shift_count <- select_samples[, .(count=.N), by=.(shift_status, comparison)]

  ## Set sample order if required.
  if (!all(samples == "all")) {
    shift_count[, comparison := factor(comparison, levels=samples)]
  }

  ## Return a table if requested.
  if (return_table) return(as.data.frame(shift_count))

  ## bar plot of shifting status.
  p <- ggplot(shift_count, aes(x=.data$comparison, y=.data$count, fill=.data$shift_status)) +
    geom_col(position="stack", ...)

  return(p)
}
