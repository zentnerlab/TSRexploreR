#' Normalize TSS counts
#'
#' @description
#' edgeR, DESeq2, or CPM normalization of TSSs or gene/transcript counts.
#'
#' @inheritParams common_params
#' @param data_type Whether TSS or gene/transcript counts should be normalized.
#' @param normalization_method Either 'edgeR', 'DESeq2', or 'CPM'.
#' @param n_samples Filter out TSSs or features not meeting the the selected threshold
#'   in this number of samples.
#'
#' @details
#' This function performs one of three normalizations on TSS or gene/transcript counts. The
#' simplest of these is counts per million (CPM), which accounts for sequencing depth. While
#' CPM is appropriate for comparing replicates, it is considered to be too simple for cases
#' in which there are expected to be substantial differences in RNA composition between samples.
#' For between-sample comparisons, the trimmed median of M-values (TMM) or median-of-ratios (MOR)
#' approaches, implemented in edgeR and DESeq2, respectively, can be used. Both of these methods
#' are designed to reduce the impact of library size on such comparisons. Prior to TMM or MOR
#' normalization, it is recommend to remove features with few or no reads, as they may bias the 
#' final results. To facilitate this filtering, two arguments are provided: 'threshold' and 
#' 'n_samples'. Features must have greater than or equal to 'threshold' number of raw counts 
#' in at least 'n_samples' number of samples to proceed through normalization.
#'
#' @return TSRexploreR object with normalized counts.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#' exp <- format_counts(exp, data_type="tss")
#' exp <- normalize_counts(exp, data_type="tss", method="DESeq2")
#'
#' @seealso \code{\link{plot_correlation}} for correlation analysis and visualization.
#'
#' @rdname normalize_counts-function
#' @export

normalize_counts <- function(
  experiment,
  data_type=c("tss", "tss_features"),
  method="DESeq2",
  threshold=1,
  n_samples=1
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tss_features")
  )
  method <- match.arg(str_to_lower(method), c("edger", "deseq2", "cpm"))
  assert_that(is.count(threshold))
  assert_that(is.count(n_samples))
  if (method == "DESeq2") {
    assert_that(!is.null(experiment@meta_data$sample_sheet))
  }

  ## Get selected samples.
  select_samples <- extract_counts(experiment, data_type, "all")

  ## Generate count matrix from samples.
  select_samples <- .count_matrix(select_samples, "tss")

  ## Filter counts.
  filtered_counts <- select_samples[
    rowSums(select_samples >= threshold) >= n_samples,
  ]

  ## If method is DESeq2, prepare coldata.
  if (method == "deseq2") {
    coldata <- experiment@meta_data$sample_sheet
    coldata <- coldata[match(coldata$sample_name, colnames(filtered_counts)), ]
  }

  ## Normalize filtered counts.
  normalized_counts <- switch(
    method,
    "edger"=.edger_normalize(filtered_counts),
    "cpm"=cpm(filtered_counts),
    "deseq2"=.deseq2_normalize(filtered_counts, coldata)
  )

  ## Add normalized counts to sample tables.
  normalized_counts <- normalized_counts %>%
    as.data.table(keep.rownames="FHASH") %>%
    melt(
      id.vars="FHASH", variable.name="sample",
      value.name="normalized_score"
    )

  sample_tables <- experiment %>%
    extract_counts(data_type, "all") %>%
    rbindlist(idcol="sample") %>%
    {merge(normalized_counts, ., by=c("sample", "FHASH"), all.y=TRUE)} %>%
    as_granges %>%
    sort %>%
    as.data.table %>%
    split(by="sample", keep.by=FALSE)

  experiment <- set_count_slot(
    experiment, sample_tables, "counts",
    data_type, "raw"
  )
  
  return(experiment)
}

#' edgeR Normalization
#'
#' @param count_matrix count matrix

.edger_normalize <- function(count_matrix) {
  
  ## Input check.
  assert_that(is.matrix(count_matrix))

  ## TMM Normalization.
  normalized_counts <- count_matrix %>%
    DGEList %>%
    calcNormFactors %>%
    cpm

  return(normalized_counts)   
}

#' DESeq2 Normalization
#'
#' @param count_matrix count matrix.
#' @param design design formula.
#' @param coldata column data.
#' @param blind Whether normalization is blind.

.deseq2_normalize <- function(count_matrix, coldata) {

  ## Validate inputs.
  assert_that(is.matrix(count_matrix))
  assert_that(is.data.frame(coldata))

  ## DESeq2 normalization.
  normalized_counts <- count_matrix %>%
    DESeqDataSetFromMatrix(colData=coldata, design= ~ 1) %>%
    estimateSizeFactors %>%
    counts(normalized=TRUE)

  return(normalized_counts)
}
