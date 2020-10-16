
#' Annotate Data
#'
#' @description
#' Use the ChIPseeker package to annotate TSSs or TSRs relative to known genes or transcripts.
#'
#' @import tibble
#' @importFrom ChIPseeker annotatePeak
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param annotation_data Path to annotation file or loaded TxDb object
#' @param data_type Whether to annotate TSSs or TSRs
#' @param feature_type Annotate at the gene or transcript level
#' @param upstream Bases upstream of TSS for 'promoter' annotation
#' @param downstream Bases downstream of TSS for 'promoter' annotation
#'
#' @details
#' This function attempts to assign TSSs or TSRs to the nearest genomic feature.
#' Genomic annotation data can be provided as either a 'GTF' or 'GFF' file,
#'   or as a TxDb package from bioconductor.
#'
#' 'feature_type' allows to you link TSSs or TSRs to genes or transcripts.
#' Furthermore, the size of the promoter region can be defined using
#'   'upstream' and 'downstream', which are relative to the TSSs
#'   defined in your annotation data.
#'
#' @return tsrexplorer object with annotated features
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data = annotation,
#'   data_type = "tss", feature_type = "transcript"
#' )
#'
#' @rdname annotate_features-function
#' @export

annotate_features <- function(
  experiment,
  annotation_data,
  data_type = c("tss", "tsr", "tss_diff", "tsr_diff"),
  feature_type = c("gene", "transcript"),
  upstream = 1000,
  downstream = 100
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is(annotation_data, "character") | is(annotation_data, "TxDb"),
    msg="annotation_data must be an annotation file or TxDb object"
  )
  data_type <- match.arg(data_type, c("tss", "tsr")) 
  feature_type <- match.arg(feature_type, c("gene", "transcript"))
  assert_that(is.count(upstream) && upstream > 0)
  assert_that(is.count(downstream) && downstream > 0)

  ## Load GTF.
  anno_type <- case_when(
    is(annotation_data, "character") ~ "character",
    is(annotation_data, "TxDb") ~ "txdb"
  )

  genome_annotation <- switch(
    anno_type,
    "character"=makeTxDbFromGFF(annotation_data),
    "txdb"=annotation_data
  )

  ## Grab data from proper slot.
  counts <- switch(data_type,
    "tss"=experiment@counts$TSSs$raw,
    "tsr"=experiment@counts$TSRs$raw,
    "tss_diff"=experiment@diff_features$TSSs$results,
    "tsr_diff"=experiment@diff_features$TSRs$results
  )

  ## Annotate features.
  
  counts_annotated <- counts %>%
    map(function(x) {
      x <- as_granges(x)

      annotated <- x %>%
        annotatePeak(
          tssRegion = c(-upstream, downstream),
          TxDb = genome_annotation,
          sameStrand = TRUE,
          level = feature_type
        ) %>%
        as.data.table
      
      annotated[,
        simple_annotations := case_when(
          annotation == "Promoter" ~ "Promoter",
          str_detect(annotation, pattern="(Exon|UTR)") ~ "Exon",
          str_detect(annotation, pattern="Intron") ~ "Intron",
          str_detect(annotation, pattern="Downstream") ~ "Downstream",
          annotation == "Distal Intergenic" ~ "Intergenic"
        )
      ]
      
      return(annotated)
    })

  ## Place annotated features back into the tsrexplorer object.

  if (data_type %in% c("tss", "tsr")) {
      experiment <- set_count_slot(experiment, counts_annotated, "counts", data_type, "raw")
  } else if (data_type %in% c("tss_diff", "tsr_diff")) {
      experiment <- set_count_slot(experiment, counts_annotated, "diff_features", data_type, "results")
  }

  ## Save annotation settings to the tsrexplorer object.
  anno_settings <- data.table(
    "feature_type" = feature_type,
    "upstream" = upstream,
    "downstream" = downstream
  )
  experiment@settings[["annotation"]] <- anno_settings

  return(experiment)
}
