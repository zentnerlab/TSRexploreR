#' @import data.table
#' @importFrom stats as.formula cor ecdf median model.matrix
#'   p.adjust rbinom

NULL

#' @importFrom tibble tibble as_tibble column_to_rownames
#' @importFrom stringr
#'   str_detect str_replace str_c str_split str_to_lower
#'   str_extract str_sub str_pad
#' @importFrom purrr
#'   imap map walk iwalk discard pmap map2 flatten
#'   map_chr map_df map_lgl
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom ggplot2
#'   labs ylab xlab ylim xlim
#'   ggplot aes coord_flip
#'   facet_grid facet_wrap
#'   scale_fill_viridis_c scale_fill_viridis_d
#'   scale_color_viridis_c scale_color_viridis_d
#'   scale_fill_continuous scale_fill_manual scale_color_manual
#'   scale_x_continuous scale_x_discrete scale_color_continuous
#'   theme_bw theme_minimal
#'   theme element_text element_blank unit margin element_rect
#'   geom_density geom_point geom_col geom_vline geom_tile geom_line
#'   geom_violin geom_boxplot geom_jitter geom_raster geom_histogram
#'   geom_bar geom_hline
#' @importFrom forcats fct_rev fct_inorder fct_reorder fct_relevel fct_recode
#' @importFrom dplyr pull case_when desc
#' @importFrom assertthat
#'   assert_that is.count is.flag is.string has_name has_attr
#'   is.readable
#' @importFrom rlang .data enquo enquos quo_is_null "!!" "!!!" quo_text ensym

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
