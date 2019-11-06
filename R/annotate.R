
#' Annotate Data
#'
#' Use ChIPseeker package to annotate TSSs or TSRs.
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom ChIPseeker annotatePeak
#' @importFrom purrr map
#' @importFrom GenomicFeatures makeTxDbFromGFF
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
		annotated <- experiment@experiment$TSSs
	} else if (data_type == "tsr") {
		annotated <- experiment@experiment$TSRs
	}
		
	## Annotate features.
	annotated <- map(annotated,
		~ annotatePeak(.,
			tssRegion = c(-upstream, downstream),
			TxDb = genome_annotation,
			sameStrand = TRUE,
			level = feature_type,
		) %>% as_tibble(.name_repair="unique")
	)

	## Place annotated features back into tsrexplorer object.
	if (data_type == "tss") {
		experiment@annotated$TSSs <- annotated
	} else if (data_type == "tsr") {
		experiment@annotated$TSRs <- annotated
	}

	return(experiment)
}
