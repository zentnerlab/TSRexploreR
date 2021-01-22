
#' Shifting Rank Plot
#'
#' @inheritParams common_params
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
#' @inheritParams common_params
#' @param ... Arguments passed to geom_col
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
