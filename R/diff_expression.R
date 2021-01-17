#' Model for differential feature analysis
#'
#' Find differential TSSs, TSRs, or genes/transcripts.
#'
#' @inheritParams common_params
#' @param data_type Whether TSS, TSR, or gene/transcript counts should be analyzed.
#' @param formula DE formula.
#' @param method Either 'DESeq2' or 'edgeR'.
#'
#' @return DGEList object with fitted model.
#'
#' @export

fit_de_model <- function(
  experiment,
  formula,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  samples="all",
  method="DESeq2"
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))
  assert_that(is.character(samples) && (all(samples == "all") || length(samples) >= 6))
  assert_that(is(formula, "formula"))
  method <- match.arg(str_to_lower(method), c("deseq2", "edger"))

  ## Design table.
  sample_sheet <- copy(experiment@meta_data$sample_sheet)
  sample_sheet[, c("file_1", "file_2") := NULL]
  assert_that(all(all.vars(formula) %in% colnames(sample_sheet)))
  sample_sheet <- column_to_rownames(sample_sheet, "sample_name")

  ## Get data from appropriate slot and convert to count matrix.
  sample_data <- experiment %>%
    extract_counts(data_type, samples) %>%
    .count_matrix(data_type)

  ## Ensure rows of sample sheet match columns of count matrix.
  sample_sheet <- sample_sheet[
    match(colnames(sample_data), rownames(sample_sheet)),
    , drop=FALSE
  ]

  ## Build the DE model.
  fitted_model <- switch(
    method,
    "edger"=.edger_model(sample_data, sample_sheet, formula),
    "deseq2"=.deseq2_model(sample_data, sample_sheet, formula)
  )

  ## Store model in TSRexploreR object.
  if (data_type == "tss") {
    experiment@diff_features$TSSs$model <- fitted_model
  } else if (data_type == "tsr") {
    experiment@diff_features$TSRs$model <- fitted_model
  } else if (data_type == "tss_features") {
    experiment@diff_features$TSS_features$model <- fitted_model
  } else if (data_type == "tsr_features") {
    experiment@diff_features$TSR_features$model <- fitted_model
  }

  return(experiment)
}

#' edgeR Differential Expression Model
#'
#' @inheritParams common_params
#' @param count_data Count matrix.
#' @param formula Differential expression formula.

.edger_model <- function(
  count_data,
  sample_sheet,
  formula
) {

  ## Check inputs.
  assert_that(is.matrix(count_data))
  assert_that(is.data.frame(sample_sheet))
  assert_that(is(formula, "formula"))

  ## Design matrix.
  design <- model.matrix(formula, data=sample_sheet)

  ## Differential Expression.
  de_model <- count_data %>%
    DGEList(samples=sample_sheet) %>%
    {.[
      filterByExpr(.,
        design,
        min.count=3,
        min.total.count=9
      ), ,
      keep.lib.sizes=FALSE
    ]} %>%
    calcNormFactors %>%
    estimateDisp(design) %>%
    glmQLFit(design)

  return(de_model)

}

#' DESeq2 Differential Expression Model
#'
#' @inheritParams common_params
#' @param count_data Count matrix.
#' @param formula Differential expression formula.

.deseq2_model <- function(
  count_data,
  sample_sheet,
  formula
) {

  ## Check inputs.
  assert_that(is.matrix(count_data))
  assert_that(is.data.frame(sample_sheet))
  assert_that(is(formula, "formula"))

  ## Differential expression.
  de_model <- count_data %>%
    DESeqDataSetFromMatrix(colData=sample_sheet, design=formula) %>%
    DESeq

  return(de_model)
}

#' Analyze Differential Features
#'
#' Find differential TSSs, TSRs, or features from edgeR or DESeq2 model.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom edgeR glmQLFTest
#' @importFrom purrr map_dbl
#'
#' @inheritParams common_params
#' @param data_type Whether the input was generated from TSSs, TSRs, or genes/transcripts
#' @param comparison_name The name given to the comparison when stored back into the TSRexploreR object.
#' @param comparison_type For DEseq2, either 'contrast' or 'name'. For edgeR, either 'contrast' or 'coef'.
#' @param comparison For DESeq2, the contrast or name. For edgeR, the coefficients or contrasts.
#' @param shrink_lfc For DESeq2, whether the log2 fold changes are shrunk (TRUE) or not (FALSE).
#' 
#' @return tibble of differential TSRs.
#'
#' @rdname differential_expression-function
#'
#' @export

differential_expression <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  comparison_name,
  comparison_type,
  comparison,
  shrink_lfc=FALSE
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))
  assert_that(is.string(comparison_name))
  comparison_type <- match.arg(
    str_to_lower(comparison_type),
    c("name", "contrast", "coef")
  )
  assert_that(is.vector(comparison))
  assert_that(is.flag(shrink_lfc))

  ## Get appropriate model.
  de_model <- switch(
    data_type,
    "tss"=experiment@diff_features$TSSs$model,
    "tsr"=experiment@diff_features$TSRs$model,
    "tss_features"=experiment@diff_features$TSS_features$model,
    "tsr_features"=experiment@diff_features$TSR_features$model
  )

  ## Retrieve the DE method.
  de_method <- case_when(
    is(de_model, "DESeqDataSet") ~ "deseq2",
    is(de_model, "DGEGLM") ~ "edger"
  )

  ## Run differential expression.
  de_args <- list()
  if (de_method == "edger") {
    de_args[[comparison_type]] <- comparison
    de_results <- do.call(glmQLFTest, c(list(de_model), de_args))
  } else if (de_method == "deseq2") {
    if (shrink_lfc) {
      de_args <- list(type="apeglm", coef=comparison)
      de_results <- do.call(lfcShrink, c(list(de_model), de_args))
    } else {
      de_args[[comparison_type]] <- comparison
      de_args[["cooksCutoff"]] <- FALSE
      de_results <- do.call(results, c(list(de_model), de_args))
    }
  }

  ## Get table of results.
  if (de_method == "deseq2") {
    de_results <- as.data.table(de_results, keep.rownames="feature")
  } else if (de_method == "edger") {
    de_results <- as.data.table(de_results$table, keep.rownames="feature")
  }

  if (de_method == "deseq2") {
    de_results[, lfcSE := NULL]
    setnames(
      de_results,
      old=c("log2FoldChange", "baseMean"),
      new=c("log2FC", "mean_expr")
    )
  } else if (de_method == "edger") {
    de_results[, F := NULL]
    setnames(
      de_results,
      old=c("logFC", "PValue", "logCPM"),
      new=c("log2FC", "pvalue", "mean_expr")
    )
    de_results[, padj := p.adjust(pvalue, method="fdr")]
  }

  ## Split out ranges.
  de_results[,
    c("seqnames", "start", "end", "strand") :=
    tstrsplit(feature, ":")
  ][,
    c("start", "end") := lapply(.SD, as.numeric),
    .SDcols=c("start", "end")
  ]

  ## Add differential expression data back to TSRexploreR object.
  if (data_type == "tss") {
    experiment@diff_features$TSSs$results[[comparison_name]] <- de_results
  } else if (data_type == "tsr") {
    experiment@diff_features$TSRs$results[[comparison_name]] <- de_results
  } else if (data_type == "tss_features") {
    experiment@diff_features$TSS_features$results[[comparison_name]] <- de_results
  } else if (data_type == "tsr_features") {
    experiment@diff_features$TSR_features[[comparison_name]] <- de_results
  }

  return(experiment)
}

#' Mark DE Status
#'
#' @inheritParams common_params
#' @param de_results Results of DE.

.de_status <- function(
  de_results,
  log2fc_cutoff,
  fdr_cutoff
) {

  ## Check inputs.
  assert_that(is.data.frame(de_results))
  assert_that(is.numeric(log2fc_cutoff) && log2fc_cutoff >= 0)
  assert_that(is.numeric(fdr_cutoff) && (fdr_cutoff > 0 & fdr_cutoff <= 1))

  ## Mark DE status.
  de_results[,
    de_status := case_when(
      is.na(padj) | is.na(log2FC) ~ "unchanged",
      padj > fdr_cutoff | abs(log2FC) < log2fc_cutoff ~ "unchanged",
      padj <= fdr_cutoff & log2FC >= log2fc_cutoff ~ "up",
      padj <= fdr_cutoff & log2FC <= -log2fc_cutoff ~ "down"
    )
  ][,
    de_status := factor(de_status, levels=c("up", "unchanged", "down"))
  ]

}

#' DE Table
#'
#' Output a table with differential features
#'
#' @inheritParams common_params
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The name of the DE comparison
#' @param de_type A single value or combination of 'up, 'unchanged', and/or 'down' (qq a list?)
#'
#' @rdname de_table-function
#' @export

de_table <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  de_comparisons="all",
  de_type=c("up", "unchanged", "down")
) {
  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))
  assert_that(is.character(de_comparisons))
  de_type <- match.arg(str_to_lower(de_type), c("up", "unchanged", "down"), several.ok=TRUE)

  ## Grab tables.
  de_tables <- experiment %>%
    extract_de(data_type, de_comparisons) %>%
    bind_rows

  ## Filter tables.
  de_tables <- de_tables[DE %in% de_type]

  ## Return tables.
  return(de_tables)
}
