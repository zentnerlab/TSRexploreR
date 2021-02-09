#' Annotate Features
#'
#' @description
#' Use the ChIPseeker package to annotate TSSs or TSRs relative to known genes or transcripts.
#'
#' @inheritParams common_params
#' @param data_type Whether to annotate TSSs ('tss') or TSRs ('tsr').
#' @param feature_type Whether to annotate at the 'gene' or 'transcript' level.
#' @param upstream Number of bases upstream of TSS for 'promoter' annotation (integer).
#' @param downstream Number of bases downstream of TSS for 'promoter' annotation (integer).
#'
#' @details
#' This function attempts to assign TSSs or TSRs to the nearest genomic feature.
#' Genomic annotation data can be provided as either a 'GTF' or 'GFF' file,
#'   or as a TxDb package from Bioconductor.
#'
#' 'feature_type' allows to you link TSSs or TSRs to genes or transcripts.
#' Furthermore, the size of the promoter region can be defined using
#'   'upstream' and 'downstream', which are relative to the TSSs
#'   defined in your annotation data.
#' TSSs or TSRs overlapping a gene on the opposite strand
#'   will be marked as 'Antisense'.
#'
#' @return TSRexploreR object with annotation data added to TSS or TSR tables.
#'
#' @examples
#' data(TSSs_reduced)
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#'
#' exp <- TSSs %>%
#'   tsr_explorer(genome_annotation=annotation) %>%
#'   format_counts(data_type="tss")
#'
#' # Annotate TSSs
#' exp <- annotate_features(exp, data_type="tss")
#'
#' @export

annotate_features <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_diff", "tsr_diff", "shift"),
  feature_type=c("gene", "transcript"),
  genome_annotation=NULL,
  upstream=1000,
  downstream=100
) {

  ## Check if ChIPseeker is installed.
  if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
    stop("Package \"ChIPseeker\" needed for this function to work. Please install it.",
      call. = FALSE)
  }

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is.null(genome_annotation) ||
    is(genome_annotation, "character") || is(genome_annotation, "TxDb"),
    msg="genome_annotation must be a GTF/GFF3 annotation file or TxDb object"
  )
  data_type <- match.arg(data_type, c("tss", "tsr", "tss_diff", "tsr_diff", "shift")) 
  feature_type <- match.arg(feature_type, c("gene", "transcript"))
  assert_that(is.count(upstream))
  assert_that(is.count(downstream))

  ## Load GTF.
  genome_annotation <- .prepare_annotation(genome_annotation, experiment)

  ## Get data from proper slot.
  counts <- switch(data_type,
    "tss"=experiment@counts$TSSs$raw,
    "tsr"=experiment@counts$TSRs$raw,
    "tss_diff"=experiment@diff_features$TSSs$results,
    "tsr_diff"=experiment@diff_features$TSRs$results,
    "shift"=experiment@shifting$results
  )

  ## Annotate features.
  
  counts_annotated <- map(counts, function(x) {
    x <- .annotate(
      x, upstream, downstream, feature_type,
      genome_annotation
    )
    return(x)
  })

  ## Place annotated features back into the TSRexploreR object.

  if (data_type %in% c("tss", "tsr")) {
    experiment <- set_count_slot(experiment, counts_annotated, "counts", data_type, "raw")
  } else if (data_type %in% c("tss_diff", "tsr_diff")) {
    experiment <- set_count_slot(experiment, counts_annotated, "diff_features", data_type, "results")
  } else if (data_type == "shift") {
    experiment@shifting$results <- counts_annotated
  }

  ## Save annotation settings to the TSRexploreR object.
  anno_settings <- data.table(
    "feature_type"=feature_type,
    "upstream"=upstream,
    "downstream"=downstream
  )
  experiment@settings$annotation <- anno_settings

  return(experiment)
}

#' Annotate Features
#'
#' @inheritParams annotate_features
#' @param sample_table Sample table
#' @param annotation_data Genome annotation.

.annotate <- function(
  sample_table,
  upstream,
  downstream,
  feature_type,
  annotation_data
) {

  ## Convert to GRanges.
  sample_table <- as_granges(sample_table)

  ## Annotate.
  annotated <- sample_table %>%
    ChIPseeker::annotatePeak(
      tssRegion=c(-upstream, downstream),
      TxDb=annotation_data,
      sameStrand=TRUE,
      level=feature_type,
      verbose=FALSE
    ) %>%
    as.data.table
  annotated[, geneStrand := ifelse(geneStrand == 1, "+", "-")]

  ## Create a column with simplified annotations.
  annotated[,
    simple_annotations := case_when(
      annotation == "Promoter" ~ "Promoter",
      str_detect(annotation, pattern="(Exon|UTR)") ~ "Exon",
      str_detect(annotation, pattern="Intron") ~ "Intron",
      str_detect(annotation, pattern="Downstream") ~ "Downstream",
      annotation == "Distal Intergenic" ~ "Intergenic"
    )
  ]

  ## Mark antisense features.
  annotated[,
    simple_annotations := ifelse(
      simple_annotations != "Intergenic" & strand != geneStrand,
      "Antisense", simple_annotations
    )
  ]

  return(annotated)

}
