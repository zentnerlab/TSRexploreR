
#' Annotate Data
#'
#' Use ChIPseeker package to annotate TSSs or TSRs.
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom ChIPseeker annotatePeak
#' @importFrom purrr map
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom tidyr separate gather
#' @importFrom dplyr filter
#'
#' @param experiment tsrexplorer object with TSS Granges
#' @param annotation_file GTF file with genomic annotations
#' @param data_type Whether to annotate TSSs or TSRs
#' @param feature_type Annotate on gene or transcript level
#' @param upstream Bases upstream of TSS
#' @param downstream Bases downstream of TSS
#' @param normalized Whether to use the normalized (TRUE) or raw (FALSE) counts
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
	downstream = 100,
	normalized = TRUE
) {
	## Load GTF.
	genome_annotation <- makeTxDbFromGFF(annotation_file, format = "gtf")

	## Grab data from proper slot.
	if (data_type == "tss" & normalized) {
		annotated <- experiment@normalized_counts$TSSs
	} else if (data_type == "tss" & !(normalized)) {
		annotated <- experiment@raw_counts$TSSs
	} else if (data_type == "tsr" & normalized) {
		annotated <- experiment@normalized_counts$TSRs
	} else {
		annotated <- experiment@raw_counts$TSRs
	}

	## Prepare data for annotation.
	annotated <- annotated %>%
		separate(position, into = c("seqnames", "start", "end", "strand"), sep = "_") %>%
		gather(-seqnames, -start, -end, -strand, key = "sample", value = "score") %>%
		filter(score > 0) %>%
		split(.$sample) %>%
		map(~ select(.x, -sample) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE))
		
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
