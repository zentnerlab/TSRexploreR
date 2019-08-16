#' TMM Normalize TSSs or TSRs
#'
#' Using edgeR to TMM normalize TSSs or TSRs
#'
#' @import tibble
#' @importFrom dplyr bind_rows mutate select group_by summarize
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom GenomicRanges GRangesList reduce findOverlaps makeGRangesFromDataFrame
#' @importFrom tidyr spread complete
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object
#' @param data_type Whether TSSs, TSRs, TSS feature counts, or RNA-seq feature counts should be normalized
#'
#' @return tibble of TMM normalized read counts
#'
#' @export
#' @rdname count_normalization-function

count_normalization <- function(experiment, data_type = c("tss", "tsr", "tss_features", "rnaseq_features")) {

	## Select proper slot of data.
	if (data_type == "tss") {
		raw_counts <- map(
			experiment@experiment$TSSs,
			~as_tibble(., .name_repair = "unique") %>%
				mutate(position = paste(seqnames, start, end, strand, sep="_")) %>%
				select(position, score)
		) %>%
			bind_rows(.id = "sample") %>%
			spread(key = sample, value = score, fill = 0)
	} else if (data_type == "tsr") {
		## Merge overlapping TSRs to get consensus
		tsr_consensus <- experiment@experiment$TSRs %>%
			as("GRangesList") %>%
			unlist %>%
			reduce(ignore.strand=FALSE) %>%
			as_tibble(.name_repair = "unique") %>%
			mutate(names = paste(seqnames, start, end, strand, sep="_")) %>%
			makeGRangesFromDataFrame(keep.extra.columns = TRUE)

		raw_counts <- map(
			names(experiment@experiment$TSRs),
			~findOverlaps(
				query = tsr_consensus,
				subject = makeGRangesFromDataFrame(experiment@experiment$TSRs[[.x]])
			)  %>%
				as_tibble(.name_repair = "unique") %>%
				mutate(
					score = experiment@experiment$TSRs[[.x]][subjectHits]$nTAGs,
					position = tsr_consensus[queryHits]$names
				) %>%
				select(-queryHits, -subjectHits) %>%
				group_by(position) %>%
				summarize(score = sum(score)) %>%
				complete(position = tsr_consensus$names, fill = list(score = 0))
		) %>%
			setNames(names(experiment@experiment$TSRs)) %>%
			bind_rows(.id = "sample") %>%
			spread(key = sample, value = score, fill = 0)
	} else if (data_type == "tss_features") {
		raw_counts <- experiment@raw_counts$TSS_features %>%
			rename(position = gene_id)
	} else if (data_type == "rnaseq_features") {
		raw_counts <- experiment@raw_counts$RNAseq_features %>%
			rename(position = gene_id)
	}

	## TMM normalize counts.
	tmm_counts <- raw_counts %>%
		column_to_rownames("position") %>%
		as.matrix %>%
		DGEList %>%
		calcNormFactors %>%
		cpm %>%
		as_tibble(.name_repair = "unique", rownames = "position")

	## Put counts into proper slots.
	if (data_type == "tss") {
		experiment@raw_counts$TSSs <- raw_counts
		experiment@normalized_counts$TSSs <- tmm_counts
	} else if (data_type == "tsr") {
		experiment@raw_counts$TSRs <- raw_counts
		experiment@normalized_counts$TSRs <- tmm_counts
	} else if (data_type == "tss_features") {
		experiment@normalized_counts$TSS_features <- tmm_counts
	} else if (data_type == "rnaseq_features") {
		experiment@normalized_counts$RNAseq_features <- tmm_counts
	}

	return(experiment)
}
