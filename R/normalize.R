#' Normalize TSS counts
#'
#' @description
#' edgeR, DESeq2, or CPM normalization of TSSs or gene/transcript counts.
#'
#' @inheritParams common_params
#' @param data_type Whether TSS or gene/transcript counts should be normalized.
#' @param method Either 'edgeR', 'DESeq2', or 'CPM'.
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
#' When clustering TSSs into TSRs using 'tss_clustering',
#'   both the raw and normalized counts will be stored in the new TSRs.
#'
#' @return TSRexploreR object with normalized counts.
#'
#' @examples
#' data(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=rep(NA, 6), file_2=rep(NA, 6),
#'   condition=c(
#'     rep("Diamide", 3),
#'     rep("Untreated", 3)
#'   )
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet) %>%
#'   format_counts(data_type="tss")
#'
#' # CPM normalization
#' tsre <- normalize_counts(tsre, method="CPM")
#'
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

  ## Check whether DESeq2 or edgeR is required.
  if (method == "deseq2") {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("Package \"DESeq2\" needed for this function to work. Please install it.",
        call. = FALSE)
    }
  } else if (method %in% c("edger", "cpm")) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package \"edgeR\" needed for this function to work. Please install it.",
        call. = FALSE)
    }
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
    "cpm"=edgeR::cpm(filtered_counts),
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
    edgeR::DGEList(.) %>%
    edgeR::calcNormFactors(.) %>%
    edgeR::cpm(.)

  return(normalized_counts)   
}

#' DESeq2 Normalization
#'
#' @param count_matrix count matrix.
#' @param coldata column data.

.deseq2_normalize <- function(count_matrix, coldata) {

  ## Validate inputs.
  assert_that(is.matrix(count_matrix))
  assert_that(is.data.frame(coldata))

  ## DESeq2 normalization.
  normalized_counts <- count_matrix %>%
    DESeq2::DESeqDataSetFromMatrix(colData=coldata, design= ~ 1) %>%
    DESeq2::estimateSizeFactors(.) %>%
    DESeq2::counts(normalized=TRUE)

  return(normalized_counts)
}
