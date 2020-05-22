
#' Annotate Data
#'
#' @description
#' Use the ChIPseeker package to annotate TSSs or TSRs relative to known genes or transcripts.
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom stringr str_detect
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
	if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsrexplorer object")

	if (!is(annotation_data, "character") & !is(annotation_data, "TxDb")) {
		stop("annotation_data must be an annotation file or TxDb object")
	}
	
        if (!is(data_type, "character") || length(data_type) > 1) {
		stop("data_type must be a character")
	}
        data_type <- str_to_lower(data_type)
        if (!data_type %in% c("tss", "tsr")) stop("data_type must be 'tss' or 'tsr'")

	if (!is(feature_type, "character") || length(feature_type) > 1) {
		stop("feature_type must be 'gene' or 'transcript'")
	}
	feature_type <- str_to_lower(feature_type)
	if (!feature_type %in% c("gene", "transcript")) {
		stop("feature_type must be 'gene' or 'transcript'")
	}

	if (!is(upstream, "numeric") | !is(downstream, "numeric")) {
		stop("upstream and downstream must be positive integers")
	}
	if (upstream %% 1 != 0 | downstream %% 1 != 0) {
		stop("upstream and downstream must be positive integers")
	}
	if (upstream < 0 | downstream < 0) stop("upstream and downstream must be positive integers")

	## Load GTF.
	if (is(annotation_data, "character")) {
		genome_annotation <- makeTxDbFromGFF(annotation_data)
	} else if (is(annotation_data, "TxDb")) {
		genome_annotation <- annotation_data
	}

	## Grab data from proper slot.
	if (data_type == "tss") {
		counts <- experiment@counts$TSSs$raw
	} else if (data_type == "tsr") {
		counts <- experiment@counts$TSRs$raw
	} else if (data_type == "tss_diff") {
		counts <- experiment@diff_features$TSSs$results
	} else if (data_type == "tsr_diff") {
		counts <- experiment@diff_features$TSRs$results
	}

	## Annotate features.
	
	counts_annotated <- counts %>%
		map(function(x) {
			x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)

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
	if (data_type == "tss") {
		experiment@counts$TSSs$raw <- counts_annotated
	} else if (data_type == "tsr") {
		experiment@counts$TSRs$raw <- counts_annotated
	} else if (data_type == "tss_diff") {
		experiment@diff_features$TSSs$results <- counts_annotated
	} else if (data_type == "tsr_diff") {
		experiment@diff_features$TSRs$results <- counts_annotated
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
