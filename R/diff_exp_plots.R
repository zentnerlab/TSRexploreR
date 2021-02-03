#' MA Plot
#'
#' @description
#' Generate an MA plot for differential TSSs, TSRs, or genes/transcripts.
#'
#' @inheritParams common_params
#' @param de_comparisons Character vector of differential expression comparisons to plot.
#' @param data_type Either 'tss' or 'tsr'.
#' @param ... Arguments passed to geom_point.
#'
#' @details
#' This function generates an MA plot of the results from differential analysis of 
#' TSSs or TSRs.
#'
#' 'de_comparisons' are the names given to the comparisons from the 'comparison_name'
#'   argument of the 'differential_expression' function.
#' 'log2fc_cutoff' and 'fdr_cutoff' are the Log2 Fold Change and FDR cutoffs used
#'   for consideration of significance in the plot.
#'
#' @return ggplot2 object of MA plot.
#'
#' @seealso
#' \code{\link{fit_de_model}} to fit a differential expression model.
#' \code{\link{differential_expression}} to find differential TSSs or TSRs.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=NA, file_2=NA,
#'   condition=c(rep("Diamide", 3), rep("Untreated", 3))
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet) %>%
#'   format_counts(data_type="tss")
#'
#' # Differential TSS MA plot.
#' diff_tss <- tsre %>%
#'   fit_de_model(data_type="tss", formula= ~condition) %>%
#'   differential_expression(
#'     exp, data_type="tss",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="name",
#'     comparison="condition_Untreated_vs_Diamide"
#'   )
#' \donttest{plot_ma(diff_tss, data_type="tss")}
#'
#' # Differential TSR MA plot.
#' diff_tsr <- tsre %>%
#'   tss_clustering(threshold=3) %>%
#'   fit_de_model(data_type="tsr", formula= ~condition) %>%
#'   differential_expression(
#'     exp, data_type="tsr",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="name",
#'     comparison="condition_Untreated_vs_Diamide"
#'   )
#' \donttest{plot_ma(diff_tsr, data_type="tsr")}
#'
#' @export

plot_ma <- function(
  experiment,
  data_type=c("tss", "tsr"),
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
  p <- ggplot(de_samples, aes(x=log2(.data$mean_expr), y=.data$log2FC, color=.data$de_status)) +
    geom_point(...) +
    geom_hline(yintercept=log2fc_cutoff, lty=2) +
    geom_hline(yintercept=-log2fc_cutoff, lty=2) +
    geom_hline(yintercept=0) +
    facet_wrap(~ sample, ncol=ncol, scales="free")

  return(p)
}

#' DE Volcano Plot
#'
#' @description
#' Generate a volcano plot for differential TSSs or TSRs.
#'
#' @inheritParams common_params
#' @param data_type either 'tss' or 'tsr'.
#' @param de_comparisons Character vector of differential expression comparisons to plot.
#' @param ... Arguments passed to geom_point.
#'
#' @details
#' This function generates a volcano plot of the results from differential analysis of 
#' TSSs or TSRs.
#'
#' 'de_comparisons' are the names given to the comparisons from the 'comparison_name'
#'   argument of the 'differential_expression' function.
#' 'log2fc_cutoff' and 'fdr_cutoff' are the Log2 Fold Change and FDR cutoffs used
#'   for consideration of significance in the plot.
#'
#' @return ggplot2 object of the volcano plot.
#'
#' @seealso
#' \code{\link{fit_de_model}} to fit a differential expression model.
#' \code{\link{differential_expression}} to find differential TSSs or TSRs.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=NA, file_2=NA,
#'   condition=c(rep("Diamide", 3), rep("Untreated", 3))
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet) %>%
#'   format_counts(data_type="tss")
#'
#' # Differential TSS volcano plot.
#' diff_tss <- tsre %>%
#'   fit_de_model(data_type="tss", formula= ~condition) %>%
#'   differential_expression(
#'     exp, data_type="tss",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="name",
#'     comparison="condition_Untreated_vs_Diamide"
#'   )
#' \donttest{plot_volcano(diff_tss, data_type="tss")}
#'
#' # Differential TSR volcano plot.
#' diff_tsr <- tsre %>%
#'   tss_clustering(threshold=3) %>%
#'   fit_de_model(data_type="tsr", formula= ~condition) %>%
#'   differential_expression(
#'     exp, data_type="tsr",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="name",
#'     comparison="condition_Untreated_vs_Diamide"
#'   )
#' \donttest{plot_volcano(diff_tsr, data_type="tsr")}
#'
#' @export

plot_volcano <- function(
  experiment,
  data_type=c("tss", "tsr"),
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
#' Export differential features for use in clusterProfiler term enrichment.
#'
#' @inheritParams common_params
#' @param data_type Whether to export genes associated with differential TSSs or
#'   TSRs.
#' @param de_comparisons Character vector of differential expression comparisons to export.
#' @param keep_unchanged Logical for inclusion of genes not significantly changed in
#'   the exported list.
#' @param anno_categories Vector of annotation categories to keep.
#'   If NULL no filtering by annotation type occurs.
#'
#' @details
#' This function outputs a data.frame that is formatted for use in the 'compateCluster'
#'   function of the clusterProfiler library.
#' Importantly, the 'geneId', 'sample', and 'de_status' columns can be used in the formula
#'   notation 'geneId ~ sample + de_status'.
#'
#' 'de_comparisons' are the names given to the comparisons from the 'comparison_name'
#'   argument of the 'differential_expression' function.
#' 'log2fc_cutoff' and 'fdr_cutoff' are the Log2 Fold Change and FDR cutoffs used
#'   for consideration of significance in the plot.
#'
#' 'keep_unchanged' controls whether genes with the category of 'unchanged'
#'   (not differentially expressed) are returned in the table also.
#' Additionally, genes can be returned based on whether they have differential
#'   features within a certain relative genomic location, such as promoter.
#' This is controlled by providing a vector annotation types to 'anno_types'
#'   which will be kept.
#'
#' @return data.frame of genes and differential expression status of TSSs or TSRs.
#'
#' @seealso
#' \code{\link{fit_de_model}} to fit a differential expression model.
#' \code{\link{differential_expression}} to find differential TSSs or TSRs.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=NA, file_2=NA,
#'   condition=c(rep("Diamide", 3), rep("Untreated", 3))
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet, genome_annotation=annotation) %>%
#'   format_counts(data_type="tss") %>%
#'   annotate_features(data_type="tss")
#'
#' # Differential TSS table.
#' diff_tss <- tsre %>%
#'   fit_de_model(data_type="tss", formula= ~condition) %>%
#'   differential_expression(
#'     exp, data_type="tss",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="name",
#'     comparison="condition_Untreated_vs_Diamide"
#'   )
#' export_for_enrichment(diff_tss, data_type="tss")
#'
#' # Differential TSR table.
#' diff_tsr <- tsre %>%
#'   tss_clustering(threshold=3) %>%
#'   annotate_features(data_type="tsr") %>%
#'   fit_de_model(data_type="tsr", formula= ~condition) %>%
#'   differential_expression(
#'     exp, data_type="tsr",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="name",
#'     comparison="condition_Untreated_vs_Diamide"
#'   )
#' export_for_enrichment(diff_tsr, data_type="tsr")
#'
#' @export

export_for_enrichment <- function(
  experiment,
  data_type=c("tss", "tsr"),
  de_comparisons="all",
  log2fc_cutoff=1,
  fdr_cutoff=0.05,
  keep_unchanged=FALSE,
  anno_categories=NULL
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr"))
  assert_that(is.character(de_comparisons))
  assert_that(is.numeric(log2fc_cutoff) && log2fc_cutoff >= 0)
  assert_that(is.numeric(fdr_cutoff) && (fdr_cutoff > 0 & fdr_cutoff <= 1))
  assert_that(is.flag(keep_unchanged))
  assert_that(is.null(anno_categories) || is.character(anno_categories))

  ## Get DE comparisons.
  de_samples <- extract_de(experiment, data_type, de_comparisons)
  de_samples <- rbindlist(de_samples, idcol="sample")

  ## Keep only selected annotation categories if provided.
  if (!is.null(anno_categories)) {
    de_samples <- de_samples[simple_annotations %in% anno_categories]
  }

  ## Mark de status.
  .de_status(de_samples, log2fc_cutoff, fdr_cutoff)

  ## Remove unchanged if requested and set factor levels.
  if (keep_unchanged) {
    de_samples[, de_status := factor(
      de_status, levels=c("up", "unchanged", "down")
    )]
  } else {
    de_samples <- de_samples[de_status != "unchanged"]
    de_samples[, de_status := factor(
      de_status, levels=c("up", "down")
    )]
  }
  
  return(as.data.frame(de_samples))
}

#' Plot DE Numbers
#'
#' @description
#' Stacked barplot of the number of differential TSSs or TSRs per comparison.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#' @param de_comparisons Character vector of differential expression comparisons to plot.
#' @param keep_unchanged Whether to include (TRUE) unchanged features in the plot.
#' @param ... Additional arguments passed to geom_col.
#'
#' @details
#' Generate a stacked barplot with the number of differential TSSs or TSRs per comparison.
#'
#' 'de_comparisons' are the names given to the comparisons from the 'comparison_name'
#'   argument of the 'differential_expression' function.
#' 'log2fc_cutoff' and 'fdr_cutoff' are the Log2 Fold Change and FDR cutoffs used
#'   for consideration of significance in the plot.
#' 'keep_unchaged' controls whether non-significant feature numbers are included in the plot.
#'
#' If 'keep_unchaged' is TRUE, a table with the numbers is returned instead of the ggplot.
#' This may be useful if exact numbers underlying the plot are required.
#'
#' @return ggplot2 object of stacked barplot.
#'   If 'return_table' is TRUE, a data.frame with differentially expressed TSS/TSR numbers
#'   are returned.
#'
#' @seealso
#' \code{\link{fit_de_model}} to fit a differential expression model.
#' \code{\link{differential_expression}} to find differential TSSs or TSRs.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=NA, file_2=NA,
#'   condition=c(rep("Diamide", 3), rep("Untreated", 3))
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet) %>%
#'   format_counts(data_type="tss")
#'
#' # Differential TSS quantities plot.
#' diff_tss <- tsre %>%
#'   fit_de_model(data_type="tss", formula= ~condition) %>%
#'   differential_expression(
#'     exp, data_type="tss",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="name",
#'     comparison="condition_Untreated_vs_Diamide"
#'   )
#' \donttest{plot_num_de(diff_tss, data_type="tss")}
#'
#' # Differential TSR quantities plot.
#' diff_tsr <- tsre %>%
#'   tss_clustering(threshold=3) %>%
#'   fit_de_model(data_type="tsr", formula= ~condition) %>%
#'   differential_expression(
#'     exp, data_type="tsr",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="name",
#'     comparison="condition_Untreated_vs_Diamide"
#'   )
#' \donttest{plot_num_de(diff_tsr, data_type="tsr")}
#'
#' @export

plot_num_de <- function(
  experiment,
  data_type=c("tss", "tsr"),
  de_comparisons="all",
  log2fc_cutoff=1,
  fdr_cutoff=0.05,
  keep_unchanged=FALSE,
  return_table=FALSE,
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
  assert_that(is.flag(keep_unchanged))
  assert_that(is.flag(return_table))

  ## Get appropriate samples.
  de_samples <- extract_de(experiment, data_type, de_comparisons)
  de_samples <- rbindlist(de_samples, idcol="samples")

  ## Mark DE status.
  .de_status(de_samples, log2fc_cutoff, fdr_cutoff)
  if (keep_unchanged) {
    de_samples[, de_status := factor(
      de_status, levels=c("up", "unchanged", "down")
    )]
  } else {
    de_samples <- de_samples[de_status != "unchanged"]
    de_samples[, de_status := factor(de_status, levels=c("up", "down"))]
  }

  ## Prepare data for plotting.
  de_samples <- de_samples[, .(count=.N), by=.(samples, de_status)]

  ## Set sample order if required.
  if (!all(de_comparisons == "all")) {
    de_samples[, samples := factor(samples, levels=de_comparisons)]
  }

  ## Return table if requested.
  if (return_table) return(as.data.frame(de_samples))

  ## Plot data.
  p <- ggplot(de_samples, aes(x=.data$samples, y=.data$count, fill=.data$de_status)) +
    geom_col(position="stack", ...) +
    theme(axis.text.x=element_text(angle=45, hjust=1))

  return(p)
}
