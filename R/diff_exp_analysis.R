#' DE MA or Volcano Plot
#'
#' @description
#' Generate MA or Volcano plot for differential TSRs or Genes (RNA-seq) - confused. I see the MA-plot code but no volcano.
#'
#' @param experiment tsrexplorer object
#' @param de_comparisons Which differential expression comparisons to plot
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param plot_type Either 'ma' or 'volcano'
#' @param ncol Number of columns for the facets
#' @param ... Arguments passed to geom_point
#'
#' @details
#' This function generates an MA plot of the results from
#'   differential analysis of TSSs, TSRs, or features.
#' It is returned as a ggplot2 object.
#'
#' @return ggplot2 object of MA plot.
#'
#' @rdname plot_de_validation-function
#' @export

plot_de_validation <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  de_comparisons="all",
  plot_type="ma",
  ncol=1,
  ...
){

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))
  assert_that(is.character(de_comparisons))
  plot_type <- match.arg(str_to_lower(plot_type), c("ma", "volcano"))
  assert_that(is.count(ncol))

  ## Get differential expression tables.
  de_samples <- extract_de(experiment, data_type, de_comparisons)

  ## Prepare DE data for plotting.
  de_samples <- rbindlist(de_samples)
  de_samples[, DE := factor(DE, levels=c("up", "unchanged", "down"))]

  ## Set sample order if required.
  if (!all(de_comparisons == "all")) {
    de_samples[, sample := factor(sample, levels=de_comparisons)]
  }

  ## MA plot of differential expression.
  if (plot_type == "ma") {
    p <- ggplot(de_samples, aes(x=.data$logCPM, y=.data$log2FC, color=.data$DE)) +
      geom_point(...) +
      theme_bw() +
      scale_color_viridis_d() +
      facet_wrap(~ sample, ncol=ncol, scales="free")
  ## Volcano plot of differential expression.
  } else if (plot_type == "volcano") {
    p <- ggplot(de_samples, aes(x=.data$log2FC, y=-log10(.data$FDR))) +
      geom_point(...) +
      theme_bw() +
      scale_color_viridis_d() +
      facet_wrap(~ sample, ncol=ncol, scales="free")
  }

  return(p)
}

#' Export to clusterProfiler
#'
#' Export DEGs for use in clusterProfiler term enrichment.
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The DE comparisons to plot
#'
#' @rdname export_for_enrichment-function
#'
#' @export

export_for_enrichment <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  de_comparisons="all" 
) {
  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features")
  )
  assert_that(is.character(de_comparisons))

  ## Get DE comparisons.
  de_data <- extract_de(experiment, data_type, de_comparisons)
  de_data <- rbindlist(de_data)

  de_data <- de_data[
    DE %in% c("up", "down"),
    .(sample, feature, log2FC, FDR, DE)
  ]
  
  return(de_data)
}

#' Plot DE Numbers
#'
#' Plot number of DE features.
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The comparisons to plot
#' @param ... Additional arguments passed to geom_col
#'
#' @rdname plot_num_de-function
#' @export

plot_num_de <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  de_comparisons="all",
  ...
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features")
  )
  assert_that(is.character(de_comparisons))

  ## Get appropriate samples.
  de_samples <- extract_de(experiment, data_type, de_comparisons)
  de_samples <- rbindlist(de_samples)

  ## prepare data for plotting.
  de_data <- de_samples[, .(count=.N), by=.(sample, DE)]
  de_data[, DE := factor(DE, levels=c("up", "unchanged", "down"))]

  ## Set sample order if required.
  if (!all(de_comparisons == "all")) {
    de_samples[, sample := factor(sample, levels=de_comparisons)]
  }

  ## Plot data.
  p <- ggplot(de_data, aes(x=.data$sample, y=.data$count, fill=.data$DE)) +
    geom_col(position="stack", ...) +
    theme_bw() +
    scale_fill_viridis_d() +
    theme(axis.text.x=element_text(angle=45, hjust=1))

  return(p)
}
