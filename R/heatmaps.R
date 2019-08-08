
#' TSS Heatmap Count Matrix
#'
#' Generate count matrix to make TSS heatmap
#'
#' @include tsrexplorer.R
#' @include annotate.R
#'
#' @import tibble
#' @importFrom dplyr select filter group_by summarize bind_rows pull mutate_if case_when
#' @importFrom magrittr %>%
#' @importFrom tidyr spread gather complete
#' 
#' @param experiment tsrexplorer object with annotated TSSs
#' @param sample Name of sample to analyze
#' @param upstream Bases upstream to consider
#' @param downstream bases downstream to consider
#' @param threshold Reads required per TSS
#' @param anno_type Whether the heatmap is built on genes or transcripts ("geneId", "transcriptId")
#'
#' @return matrix of counts for each gene/transcript and position
#'
#' @export
#' @rdname tss_count_matrix-function

tss_count_matrix <- function(
	experiment, sample, upstream = 1000,
	downstream = 1000, threshold = 1,
	anno_type = c("transcriptId", "geneId")
) {
	## Prepare data to be made into count matrix
	annotated_tss <- experiment@annotated$TSSs[[sample]] %>%
		select(anno_type, distanceToTSS, score) %>%
		filter(
			distanceToTSS >= -upstream & distanceToTSS <= downstream,
			score >= threshold
		)
	
	if (anno_type == "geneId") {
		annotated_tss <- annotated_tss %>%
			complete(geneId, distanceToTSS = -upstream:downstream, fill = list(score = 0)) %>%
			spread(key = distanceToTSS, value = score) %>%
			mutate_if(is.numeric, ~log2(. + 1))
	} else {
                annotated_tss <- annotated_tss %>%
                        complete(transcriptId, distanceToTSS = -upstream:downstream, fill = list(score = 0)) %>%
                        spread(key = distanceToTSS, value = score) %>%
                        mutate_if(is.numeric, ~log2(. + 1))
	}

	return(annotated_tss)
}
