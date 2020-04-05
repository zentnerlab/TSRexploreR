
#' Annotate Data
#'
#' Use ChIPseeker package to annotate TSSs or TSRs.
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom ChIPseeker annotatePeak
#' @importFrom purrr map
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom tidyr separate gather
#' @importFrom dplyr filter select
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges "rowRanges<-"
#' @importFrom stringr str_detect
#'
#' @param experiment tsrexplorer object with TSS Granges
#' @param annotation_data Either path and file name of annotation file, or TxDb object of annotation
#' @param data_type Whether to annotate TSSs or TSRs
#' @param feature_type Annotate on gene or transcript level
#' @param upstream Bases upstream of TSS
#' @param downstream Bases downstream of TSS
#'
#' @return Annotated TSS slot in tsrchitect object
#'
#' @rdname annotate_features-function
#'
#' @export

annotate_features <- function(
	experiment,
	annotation_data,
	data_type = c("tss", "tsr"),
	feature_type = c("gene", "transcript"),
	upstream = 1000,
	downstream = 100
) {
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
	} else if (data_type == "diff_tss") {
		counts <- experiment@diff_ex
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

	## Place annotated features back into tsrexplorer object.
	if (data_type == "tss") {
		experiment@counts$TSSs$raw <- counts_annotated
	} else if (data_type == "tsr") {
		experiment@counts$TSRs$raw <- counts_annotated
	}

	## Save annotation settings to tsrexplorer object.
	anno_settings <- data.table(
		"feature_type" = feature_type,
		"upstream" = upstream,
		"downstream" = downstream
	)
	experiment@settings[["annotation"]] <- anno_settings

	return(experiment)
}
