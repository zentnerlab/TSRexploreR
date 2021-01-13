#' Density Plots
#'
#' @description
#' Generate density plots of TSS or TSR positions relative to annotated TSSs.
#'
#' @inheritParams common_params
#' @param data_type Whether to plot TSS or TSR density.
#' @param consider_score Whether the score of each TSS or TSR score be considered
#'   in addition to its unique location.
#' @param upstream Bases upstream of plot center.
#' @param downstream Bases downstream of plot center.
#' @param color Either 'default' or a valid color format to set plot color.
#' @param ... Arguments passed to geom_density.
#'
#' @details
#' This plotting function generates a density plot of TSS or TSR positions relative 
#' to annotated TSSs. The plot is returned as a ggplot2 object.
#'
#' By default, only the TSS or TSR position is considered, effectively giving every
#' TSS or TSR a score of 1. If 'consider_score' is set to TRUE, the score of each 
#' TSS or TSR will be considered when making the plot, giving more weight to 
#' stronger TSSs or TSRs.
#'
#' The region around the annotated TSS used for plotting is controlled by 'upstream' 
#' and 'downstream', which should be positive integers.
#'
#' A set of functions to control data structure for plotting are included. 'use_normalized' 
#' will use  normalized scores, which only matters if 'consider_score' is TRUE.
#' 'threshold' defines the minimum number of raw counts a TSS or TSR must have to be 
#' considered. dominant' specifies whether only the dominant TSS or TSR (determined
#' using the 'mark_dominant' function) is considered. For TSSs, this can be either 
#' dominant TSS per TSR or gene/transcript, and for TSRs it is the dominant TSR 
#' per gene/transcript. 'data_conditionals' can be used to filter, quantile, order, 
#' and/or group data for plotting.
#'   
#' @return ggplot2 object of density plot.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#' exp <- format_counts(exp, data_type="tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#' exp <- annotate_features(
#'   exp, annotation_data=annotation,
#'   data_type="tss", feature_type="transcript"
#' )
#' plot_density(exp, data_type="tss")
#'
#' @seealso
#' \code{\link{annotate_features}} to annotate TSSs or TSRs.
#'   \code{\link{mark_dominant}} to identify dominant TSSs or TSRs.
#'
#' @rdname plot_density-function
#' @export

plot_density <- function(
  experiment,
  data_type=c("tss", "tsr"),
  samples="all",
  consider_score=FALSE,
  upstream=1000,
  downstream=1000,
  threshold=NULL,
  ncol=1,
  use_normalized=FALSE,
  dominant=FALSE,
  exclude_antisense=TRUE,
  data_conditions=NULL,
  color="default",
  ...
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr"))
  assert_that(is.character(samples))
  assert_that(is.flag(consider_score))
  assert_that(is.count(upstream))
  assert_that(is.count(downstream))
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.count(ncol))
  assert_that(is.flag(use_normalized))
  assert_that(is.flag(dominant))
  assert_that(is.null(data_conditions) || is.list(data_conditions))
  assert_that(is.flag(exclude_antisense))

  ## Assign color type.
  color_type <- case_when(
    color == "default" & data_type == "tss" ~ "#431352",
    color == "default" & data_type == "tsr" ~ "#34698c",
    TRUE ~ color
  )

  ## Pull data out of appropriate slot.
  sample_data <- experiment %>%
    extract_counts(data_type, samples, use_normalized) %>%
    preliminary_filter(dominant, threshold)

  sample_data <- map(sample_data, ~ .x[dplyr::between(distanceToTSS, -upstream, downstream)])

  ## Remove antisense TSSs/TSRs if requested.
  sample_data <- map(sample_data, ~ .x[simple_annotations != "Antisense"])

  ## Condition data.
  sample_data <- condition_data(sample_data, data_conditions)

  ## Update data if score is to be considered in addition to unique position.
  sample_data <- rbindlist(sample_data, idcol="sample")
  if (consider_score) sample_data <- sample_data[rep(seq_len(.N), score)]

  ## Set sample order if required.
  if (!all(samples == "all")) {
    sample_data[, samples := factor(samples, levels=samples)]
  }

  ## Plot densities.
  p <- ggplot(sample_data, aes(.data$distanceToTSS)) +
    geom_density(fill=color_type, color=color_type, ...) +
    labs(
      x="Position Relative to Annotated TSS",
      y="Density"
    ) +
    theme_bw()

  if (!is.null(data_conditions$grouping)) {
    p <- p + facet_grid(row_groups ~ sample)
  } else if (!is.null(data_conditions$quantiling)) {
    p <- p + facet_grid(row_quantile ~ sample)
  } else {
    p <- p + facet_wrap(~ sample, ncol=ncol)
  }

  return(p)
}
