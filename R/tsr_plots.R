#' Plot TSR Metrics
#'
#' Plot selected TSR metrics.
#'
#' @inheritParams common_params
#' @param tsr_metrics Names of metrics in TSRexploreR TSR GRanges to plot.
#' @param plot_type Output either a 'violin', 'jitter', 'box', or 'boxjitter' plot.
#' @param log2_transform Whether the metric should be log2 + 1 transformed prior to plotting.
#' @param ... Arguments passed to ggplot2 plotting functions.
#'
#' @details
#'
#' Plot any TSR metrics contained within the counts data.table.
#' Metrics can be supplied as a character vector to 'tsr_metrics',
#'   and will be optionally Log2 transformed if 'log2_transform' is TRUE.
#' Valid plot types that can be supplied to 'plot_type' include
#'   'violin', 'box', 'jitter', and 'boxjitter' (a combination of boxplot and jitterplot).
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
#' @return ggplot2 object with TSR matrix plotted.
#'
#' @seealso
#' \code{\link{tsr_metrics}} to calculate additional TSR metrics.
#'
#' @examples
#' library("magrittr")
#' data(TSSs)
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3)
#' \donttest{plot_tsr_metric(tsre, "width")}
#'
#' @export

plot_tsr_metric <- function(
  experiment,
  tsr_metrics,
  plot_type="violin",
  samples="all",
  log2_transform=FALSE,
  ncol=1,
  use_normalized=FALSE,
  dominant=FALSE,
  threshold=NULL,
  data_conditions=NULL,
  ...
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(tsr_metrics))
  plot_type <- match.arg(
    str_to_lower(plot_type),
    c("violin", "jitter", "box", "boxjitter")
  )
  assert_that(is.character(samples))
  assert_that(is.flag(log2_transform))
  assert_that(is.count(ncol))
  assert_that(is.flag(use_normalized))
  assert_that(is.flag(dominant))
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.null(data_conditions) || is.list(data_conditions))

  ## Get data.
  selected_data <- experiment %>%
    extract_counts("tsr", samples, use_normalized) %>%
    preliminary_filter(dominant, threshold)

  ## Condition the data.
  selected_data <- condition_data(selected_data, data_conditions)

  grouping_status <- case_when(
    !is.null(data_conditions$quantiling) ~ "row_quantile",
    !is.null(data_conditions$grouping) ~ "row_groups",
    TRUE ~ "none"
  )

  ## Combine data into one data table.
  selected_data <- rbindlist(selected_data, idcol="sample")
  
  ## Log2 + 1 transform data if required.
  if (log2_transform) {
    selected_data[,
      (tsr_metrics) := lapply(.SD, function(x) log2(x + 1)),
      .SDcols=tsr_metrics
    ]
    setnames(
      selected_data, old=tsr_metrics,
      new=str_c("Log2(", tsr_metrics, " + 1)")
    )
    selected_data <- as.data.table(selected_data)
  }

  ## Prepare data for plotting.
  if (log2_transform) {
    selected_data <- melt(
      selected_data,
      measure.vars=str_c("Log2(", tsr_metrics, " + 1)"),
      variable.name="metric"
    )
  } else {
    selected_data <- melt(
      selected_data,
      measure.vars=tsr_metrics,
      variable.name="metric"
    )
  }

  ## Order samples if required.
  if (!all(samples == "all")) {
    selected_data[, sample := factor(sample, levels=samples)]
  }

  ## Make plot of selected TSR metric(s).
  p <- ggplot(selected_data, aes(x=.data$sample, y=.data$value))

  if (plot_type == "violin") {
    p <- p + geom_violin(aes(fill=.data$sample), ...)
  } else if (plot_type == "box") {
    p <- p + geom_boxplot(fill=NA, aes(color=.data$sample), ...)
  } else if (plot_type == "jitter") {
    p <- p + geom_jitter(aes(color=.data$sample), ...)
  } else if (plot_type == "boxjitter") {
    p <- p +
      geom_jitter(color="lightgrey", ...) +
      geom_boxplot(fill=NA, aes(color=.data$sample), outlier.shape=NA)
  }
    
  p <- p +
    theme_bw() +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )

  if (grouping_status != "none") {
    setnames(selected_data, old=grouping_status, new="grouping")
    p <- p + facet_wrap(grouping ~ metric, ncol=ncol, scales="free")
  } else {
    p <- p + facet_wrap(~ metric, ncol=ncol, scales="free")
  }

  return(p)
}
