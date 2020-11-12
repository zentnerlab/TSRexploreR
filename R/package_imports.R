#' @import data.table

NULL

#' @importFrom tibble tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom stringr
#'   str_detect str_replace str_c str_split str_to_lower
#'   str_extract str_sub str_pad
#' @importFrom purrr imap map walk iwalk discard pmap
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom ggplot2
#'   labs ylab xlab ylim xlim
#'   ggplot aes coord_flip
#'   facet_grid facet_wrap
#'   scale_fill_viridis_c scale_fill_viridis_d
#'   scale_color_viridis_c scale_color_viridis_d
#'   scale_fill_continuous scale_fill_manual scale_color_manual
#'   scale_x_continuous scale_x_discrete
#'   theme_bw theme_minimal
#'   theme element_text element_blank unit margin element_rect
#'   geom_density geom_point geom_col geom_vline geom_tile geom_line
#'   geom_violin geom_boxplot geom_jitter geom_raster geom_histogram
#'   geom_bar
#' @importFrom rlang .data
#' @importFrom forcats fct_rev fct_inorder fct_reorder fct_relevel
#' @importFrom dplyr pull case_when
#' @importFrom assertthat assert_that is.count is.flag is.string has_name
#' @importFrom knitr kable

NULL

#' @importFrom rtracklayer import export
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom plyranges
#'   as_granges stretch join_overlap_left_directed
#'   anchor_3p anchor_5p
#'   bind_ranges reduce_ranges_directed
#' @importFrom Biostrings DNAStringSet getSeq

NULL

#' @importFrom BiocGenerics width start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom SummarizedExperiment rowRanges SummarizedExperiment rowData
#'   assay "assay<-"
#' @importFrom S4Vectors "metadata<-" metadata DataFrame mcols

NULL

#' @importFrom edgeR DGEList calcNormFactors cpm filterByExpr
#' @importFrom DESeq2
#'   DESeqDataSetFromMatrix estimateSizeFactors
#'   DESeq lfcShrink results rlog counts

NULL
