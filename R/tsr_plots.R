#' Plot TSR Metrics
#'
#' Plot selected TSR metrics.
#'
#' @param experiment TSRexploreR object with TSR GRanges
#' @param tsr_metrics Names of metrics in TSRexploreR TSR GRanges to plot
#' @param plot_type Output either a 'violin', 'jitter', 'box', or 'boxjitter' plot (default: violin)
#' @param samples Either 'all' or a vector of sample names to analyze
#' @param log2_transform Whether the metric should be log2 + 1 transformed prior to plotting
#' @param ncol Number of columns to use when plotting multiple samples
#' @param use_normalized Whether to use the CPM-normalized counts
#' @param dominant Whether to only consider dominant TSRs
#' @param threshold Keep only TSRs with at least this number of raw counts
#' @param data_conditions Condition the data (filter and quantile/group available)
#' @param ... Arguments passed to ggplot2 plotting functions
#'
#' @return ggplot2 object with TSR matrix plotted
#'
#' @rdname plot_tsr_metric-function
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
  data_conditions=NA,
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
  if (all(!is.na(data_conditions)) && !is(data_conditions, "list")) stop("data_conditions must in list form")

  ## Get data.
  selected_data <- experiment %>%
    extract_counts("tsr", samples, use_normalized) %>%
    preliminary_filter(dominant, threshold)

  ## Condition the data.
  if (all(!is.na(data_conditions))) {
    selected_data <- do.call(group_data, c(list(signal_data=selected_data), data_conditions))
  }

  groupings <- !is.na(data_conditions) & any(names(data_conditions) %in% c("grouping", "quantile_by"))

  ## Combine data into one data table.
  selected_data <- rbindlist(selected_data, idcol="sample")
  
  if (samples == "all") {
    selected_data[, sample := fct_inorder(factor(sample))]
  } else {
    selected_data[, sample := factor(sample, levels=samples)]
  }

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

  if (groupings) {
    p <- p + facet_wrap(grouping ~ metric, ncol=ncol, scales="free")
  } else {
    p <- p + facet_wrap(~ metric, ncol=ncol, scales="free")
  }

  return(p)
}
