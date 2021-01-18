#' @import data.table
#' @import ggplot2
#' @import tibble 

NULL

#' @importFrom tibble tibble as_tibble
#' @importFrom stringr
#'   str_detect str_replace str_c str_split str_to_lower
#'   str_extract str_sub str_pad
#' @importFrom purrr
#'   imap map walk iwalk discard pmap map2 flatten
#'   map_chr
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom rlang .data
#' @importFrom forcats fct_rev fct_inorder fct_reorder fct_relevel
#' @importFrom dplyr pull case_when
#' @importFrom assertthat
#'   assert_that is.count is.flag is.string has_name has_attr
#'   is.readable
#' @importFrom knitr kable

NULL

#' @importFrom rtracklayer import export
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom plyranges
#'   as_granges stretch join_overlap_left_directed mutate
#'   anchor_3p anchor_5p anchor_center
#'   bind_ranges reduce_ranges_directed
#' @importFrom Biostrings DNAStringSet getSeq writeXStringSet
#' @importFrom Rsamtools indexFa FaFile ScanBamParam scanBamFlag

NULL

#' @importFrom BiocGenerics width start end "strand<-"
#' @importFrom GenomeInfoDb seqnames seqlevels "seqlengths<-" seqlengths
#' @importFrom SummarizedExperiment rowRanges SummarizedExperiment rowData
#'   assay "assay<-"
#' @importFrom S4Vectors "metadata<-" metadata DataFrame mcols elementMetadata

NULL

#' @importFrom edgeR
#'   DGEList calcNormFactors cpm filterByExpr
#'   estimateDisp glmQLFit
#' @importFrom DESeq2
#'   DESeqDataSetFromMatrix estimateSizeFactors
#'   DESeq lfcShrink results rlog counts
#' @importFrom apeglm apeglm

NULL
