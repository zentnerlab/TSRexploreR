#' CPM Normalize Counts
#'
#' @description
#' CPM normalize the TSS, TSR, and/or feature counts.
#'
#' @param experiment TSRexploreR object
#' @param data_type 'tss', 'tsr', 'tss_features', or 'tsr_features'
#'
#' @details
#' Counts Per Million (CPM) is an approach for read number normalization, 
#'   allowing comparison between samples for various downstream analyses and plots.
#' For more quantitative comparisons, TMM normalization is recommended (this should be fleshed out more qq)
#'
#' @return TSRexploreR object with CPM values.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' tsre_exp <- cpm_normalize(tsre_exp, data_type="tss")
#'
#' @rdname cpm_normalize-function
#' @export

cpm_normalize <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features")
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))

  ## Get selected samples.
  select_samples <- switch(
    data_type,
    "tss"=experiment@counts$TSSs$raw,
    "tsr"=experiment@counts$TSRs$raw,
    "tss_features"=experiment@counts$TSS_features$raw,
    "tsr_features"=experiment@counts$TSR_features$raw
  )

  ## CPM normalize counts.
  cpm_counts <- select_samples %>%
    map(function(x) {
      x[, cpm := cpm(score)]
      return(x)
    })

  ## Add CPM-normalized counts back to TSRexploreR object.
  experiment <- set_count_slot(
    experiment, cpm_counts,
    "counts", data_type, "raw"
  )

  return(experiment)
}

#' TMM Normalize TSSs or TSRs
#'
#' @description
#' Using edgeR to TMM normalize TSSs or TSRs.
#'
#' @param experiment TSRexploreR object
#' @param data_type Whether TSSs, TSRs, or RNA-seq & 5' feature counts should be normalized
#' @param normalization_method Either 'edgeR', 'DESeq2', or 'CPM'
#' @param threshold Consider only features with at least this number of raw counts
#' @param n_samples Filter out positions without features meeting the the selected threshold
#'   in this number of samples
#'
#' @details
#' The TMM normalization method, employed by edgeR, is designed to reduce the influence of
#'   library size on inter-sample comparisons.
#'
#' For TMM normalization, it is recommended to remove features with few or no reads as
#'   these may bias the final results.
#' To facilitate this filtering, two arguments are provided, 'threshold' and 'n_samples'.
#' Features must have greater than or equal to 'threshold' number of raw counts in at least
#'   'n_samples' number of samples to be retained.
#'
#' @return TSRexploreR object with tmm normalized count matrices
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' tsre_exp <- count_matrix(tsre_exp, data_type="tss")
#' tsre_exp <- tmm_normalize(exp, data_type="tss")
#'
#' @seealso \code{\link{count_matrix}} to prepare the matrices.
#'   \code{\link{plot_correlation}} for various correlation plots.
#'
#' @rdname normalize_counts-function
#' @export

normalize_counts <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  method="DESeq2",
  threshold=1,
  n_samples=1
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features")
  )
  method <- match.arg(str_to_lower(method), c("edger", "deseq2", "cpm"))
  assert_that(is.count(threshold))
  assert_that(is.count(n_samples))

  ## Get selected samples.
  select_samples <- extract_matrix(experiment, data_type, "all")

  ## Filter counts.
  sample_matrix <- as.data.table(assay(select_samples, "counts"))
  sample_matrix[, match := rowSums((.SD >= threshold)) >= n_samples]
  keep_ids <- which(sample_matrix[, match])
  
  select_samples <- select_samples[keep_ids, ]
  filtered_counts <- assay(select_samples, "counts")

  ## If method is DESeq2, prepare coldata.
  if (method == "deseq2") {
    coldata <- experiment@meta_data$sample_sheet
    coldata <- switch(
      data_type,
      "tss"=coldata[match(coldata$tss_name, colnames(filtered_counts)), ],
      "tsr"=coldata[match(coldata$tsr_name, colnames(filtered_counts)), ]
    )
  }

  ## Normalize filtered counts.
  normalized_counts <- switch(
    method,
    "edger"=.edger_normalize(filtered_counts),
    "cpm"=cpm(filtered_counts),
    "deseq2"=.deseq2_normalize(filtered_counts, coldata)
  )

  ## Create filtered and normalized RangedSummarizedExperiment.
  assay(select_samples, "normalized") <- normalized_counts

  ## Add normalized count matrix to TSRexploreR object.
  experiment <- set_count_slot(
    experiment, select_samples,
    "counts", data_type, "matrix"
  )

  ## Add normalized counts to sample tables.
  if (data_type != "tsr") {
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
  }
  
  return(experiment)
}

#' edgeR Normalize
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

#' DESeq2 Normalize
#'
#' @param count_matrix count matrix
#' @param design design formula
#' @param coldata column data
#' @param blind Whether normalization is blind

.deseq2_normalize <- function(count_matrix, coldata) {

  ## Validate inputs.
  assert_that(is.matrix(count_matrix))
  assert_that(is.data.frame(coldata))

  ## DESeq2 Normalization.
  normalized_counts <- count_matrix %>%
    DESeqDataSetFromMatrix(colData=coldata, design= ~ 1) %>%
    estimateSizeFactors %>%
    counts(normalized=TRUE)

  return(normalized_counts)
}
