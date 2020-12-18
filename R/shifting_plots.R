
#' Shifting Rank Plot
#'
#' @param experiment tsrexplorer object
#' @param samples Shifting results to plot
#' @param score_order Either descending or ascending
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
    samples <- experiment@shifting$results
  } else {
    samples <- experiment@shifting$results[samples]
  }

  samples <- rbindlist(samples, idcol="comparison")

  ## Order features by selected order.
  samples[, FID := as.character(seq_len(nrow(samples)))]

  if (score_order == "descending") {
    samples[, FID := fct_reorder(FID, shift_score, .desc=TRUE)]
  } else {
    samples[, FID := fct_reorder(FID, shift_score, .desc=FALSE)]
  }

  ## Generate the plot.
  p <- ggplot(samples, aes(x=FID, y=shift_score, fill=shift_score)) +
    geom_col() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    ) +
    facet_wrap(~ comparison, ncol=ncol, scales="free")

  return(p)
}

#' Shift Count Plot
#'
#" @param experiment tsr explorer object
#' @param samples Samples to plot
#' @param ncol Number of columns
#'
#' @export

plot_shift_count <- function(
  experiment,
  samples="all",
  ncol=3
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(is.count(ncol))

  ## Get samples.
  if (samples == "all") {
    samples <- experiment@shifting$results
  } else {
    samples <- experiment@shifting$results[samples]
  }

  samples <- rbindlist(samples, idcol="comparison")

  ## Annotate shifting status.
  samples[, shift_status := case_when(
    shift_score < 0 ~ "upstream",
    shift_score > 0 ~ "downstream",
    TRUE ~ "n.s."
  )]

  ## Get number of shifts per shifting status.
  shift_count <- samples[, .(count=.N), by=.(shift_status, comparison)]

  ## bar plot of shifting status.
  p <- ggplot(shift_count, aes(x=comparison, y=count, fill=shift_status)) +
    geom_col(position="stack")

  return(p)
}
