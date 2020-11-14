#' DE MA Plot
#'
#' @description
#' Generate a MA plot for differential TSRs or Genes (RNA-seq) - confused. I see the MA-plot code but no volcano.
#'
#' @param experiment TSRexploreR object
#' @param de_comparisons Which differential expression comparisons to plot
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param ncol Number of columns for the facets
#' @param log2fc_cutoff Log2FC cutoff value
#' @param fdr_cutoff FDR cutoff value
#' @param ... Arguments passed to geom_point
#'
#' @details
#' This function generates an MA plot of the results from
#'   differential analysis of TSSs, TSRs, or features.
#' It is returned as a ggplot2 object.
#'
#' @return ggplot2 object of MA plot.
#'
#' @export

plot_ma <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  de_comparisons="all",
  ncol=1,
  log2fc_cutoff=1,
  fdr_cutoff=0.05,
  ...
){

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features")
  )
  assert_that(is.character(de_comparisons))
  assert_that(is.count(ncol))
  assert_that(is.numeric(log2fc_cutoff) && log2fc_cutoff >= 0)
  assert_that(is.numeric(fdr_cutoff) && (fdr_cutoff <= 1 & fdr_cutoff > 0))

  ## Get differential expression tables.
  de_samples <- extract_de(experiment, data_type, de_comparisons)

  ## Prepare DE data for plotting.
  de_samples <- rbindlist(de_samples, idcol="sample")
  .de_status(de_samples, log2fc_cutoff, fdr_cutoff)

  ## Set sample order if required.
  if (!all(de_comparisons == "all")) {
    de_samples[, sample := factor(sample, levels=de_comparisons)]
  }

  ## MA plot of differential expression.
  p <- ggplot(de_samples, aes(x=.data$mean_expr, y=.data$log2FC, color=.data$de_status)) +
    geom_point(...) +
    geom_hline(yintercept=log2fc_cutoff, lty=2) +
    geom_hline(yintercept=-log2fc_cutoff, lty=2) +
    geom_hline(yintercept=0) +
    facet_wrap(~ sample, ncol=ncol, scales="free")

  return(p)
}

#' DE Volcano Plot
#'
#' @param experiment tsr explorer object
#' @param data_type either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The DE comparisons to plot
#' @param ncol Number of columns for plot
#' @param log2fc_cutoff Log2FC cutoff value
#' @param fdr_cutoff FDR cutoff value
#' @param ... Arguments passed to geom_point
#'
#' @export

plot_volcano <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  de_comparisons="all",
  log2fc_cutoff=1,
  fdr_cutoff=0.05,
  ncol=1,
  ...
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features")
  )
  assert_that(is.character(de_comparisons))
  assert_that(is.numeric(log2fc_cutoff) && log2fc_cutoff >= 0)
  assert_that(is.numeric(fdr_cutoff) && (fdr_cutoff < 1 & fdr_cutoff >= 0))
  assert_that(is.count(ncol))

  ## Get DE data.
  de_samples <- extract_de(experiment, data_type, de_comparisons)

  ## Prepare DE data for plotting.
  de_samples <- rbindlist(de_samples, idcol="sample")
  .de_status(de_samples, log2fc_cutoff, fdr_cutoff)
  de_samples[, c("de_status", "padj") := list(
    factor(de_status, levels=c("up", "unchanged", "down")),
    -log10(padj)
  )]

  ## Set sample order if required.
  if (!all(de_comparisons == "all")) {
    de_samples[, sample := factor(sample, levels=de_comparisons)]
  }

  ## Volcano plot.
  p <- ggplot(de_samples, aes(x=.data$log2FC, y=.data$padj)) +
    geom_point(aes(color=.data$de_status), ...) +
    geom_vline(xintercept=log2fc_cutoff, lty=2) +
    geom_vline(xintercept=-log2fc_cutoff, lty=2) +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=-log10(fdr_cutoff), lty=2) +
    facet_wrap(~ sample, ncol=ncol, scales="free")

  return(p) 

}

#' Export to clusterProfiler
#'
#' Export DEGs for use in clusterProfiler term enrichment.
#'
#' @param experiment TSRexploreR object
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
#' @param experiment TSRexploreR object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The comparisons to plot
#' @param logfc_cutoff Log2 fold change cutoff
#' @param fdr_cutoff FDR cutoff value
#' @param ... Additional arguments passed to geom_col
#'
#' @rdname plot_num_de-function
#' @export

plot_num_de <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  de_comparisons="all",
  log2fc_cutoff=1,
  fdr_cutoff=0.05,
  ...
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features")
  )
  assert_that(is.character(de_comparisons))
  assert_that(is.numeric(log2fc_cutoff) && log2fc_cutoff >= 0)
  assert_that(is.numeric(fdr_cutoff) && (fdr_cutoff > 0 & fdr_cutoff <= 1))

  ## Get appropriate samples.
  de_samples <- extract_de(experiment, data_type, de_comparisons)
  de_samples <- rbindlist(de_samples, idcol="samples")

  ## Mark DE status.
  .de_status(de_samples, log2fc_cutoff, fdr_cutoff)
  de_samples[, de_status := factor(
    de_status, levels=c("up", "unchanged", "down")
  )]

  ## prepare data for plotting.
  de_samples <- de_samples[, .(count=.N), by=.(samples, de_status)]

  ## Set sample order if required.
  if (!all(de_comparisons == "all")) {
    de_samples[, samples := factor(samples, levels=de_comparisons)]
  }

  ## Plot data.
  p <- ggplot(de_samples, aes(x=.data$samples, y=.data$count, fill=.data$de_status)) +
    geom_col(position="stack", ...) +
    theme(axis.text.x=element_text(angle=45, hjust=1))

  return(p)
}
