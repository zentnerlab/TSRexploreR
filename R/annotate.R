
#' Annotate Data
#'
#' Use ChIPseeker package to annotate TSSs or TSRs.
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom ChIPseeker annotatePeak
#' @importFrom purrr map
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr separate gather
#' @importFrom dplyr filter select
#'
#' @param experiment tsrexplorer object with TSS Granges
#' @param annotation_file GTF file with genomic annotations
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
	annotation_file,
	data_type = c("tss", "tsr"),
	feature_type = c("gene", "transcript"),
	upstream = 1000,
	downstream = 100
) {
	## Load GTF.
	genome_annotation <- makeTxDbFromGFF(annotation_file, format = "gtf")

	## Grab data from proper slot.
	if (data_type == "tss") {
		raw <- experiment@counts$TSSs$raw
		cpm_counts <- experiment@counts$TSSs$cpm
	} else if (data_type == "tsr") {
		raw <- experiment@counts$TSRs$raw
		cpm_counts <- experiment@counts$TSRs$cpm
	}

	## Annotate features.
	raw_annotated <- raw %>%
		map(
			~ annotatePeak(.,
				tssRegion = c(-upstream, downstream),
				TxDb = genome_annotation,
				sameStrand = TRUE,
				level = feature_type,
			) %>% as_tibble(.name_repair="unique")
	)

	cpm_annotated <- cpm_counts %>%
                map(
                        ~ annotatePeak(.,
                                tssRegion = c(-upstream, downstream),
                                TxDb = genome_annotation,
                                sameStrand = TRUE,
                                level = feature_type,
                        ) %>% as_tibble(.name_repair="unique")
		)


	## Place annotated features back into tsrexplorer object.
	if (data_type == "tss") {
		experiment@annotated$TSSs <- list("raw" = raw_annotated, "cpm" = cpm_annotated)
	} else if (data_type == "tsr") {
		experiment@annotated$TSRs <- list("raw" = raw_annotated, "cpm" = cpm_annotated)
	}

	return(experiment)
}
