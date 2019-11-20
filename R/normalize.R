#' TMM Normalize TSSs or TSRs
#'
#' Using edgeR to TMM normalize TSSs or TSRs
#'
#' @import tibble
#' @importFrom dplyr mutate select bind_rows group_by summarize mutate_if left_join
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom GenomicRanges GRangesList reduce findOverlaps makeGRangesFromDataFrame score mcols
#' @importFrom tidyr spread complete
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object
#' @param data_type Whether TSSs, TSRs, RNA-seq & five-prime feature counts should be normalized
#' @param threshold Filter out positions missing at least 'n_samples' number of samples with reads greater than or equal to threshold
#' @param n_samples Filter out positions missing at least n_samples number of samples with reads greater than or equal to 'threshold'
#' @param samples Vector with names of samples to include in he normalization
#'
#' @return tibble of TMM normalized read counts
#'
#' @rdname count_normalization-function
#'
#' @export

count_normalization <- function(
	experiment,
	data_type = c("tss", "tsr", "features"),
	samples = "all",
	threshold = 1,
	n_samples = 1
) {

	## Select proper slot of data.
	if (data_type == "tss") {
		if (samples == "all") samples <- names(experiment@experiment$TSSs)
		select_samples <- experiment@experiment$TSSs[samples]

		raw <- select_samples %>%
			map(function(x) {
				if ("nTAGs" %in% names(mcols(x))) score(x) <- x$nTAGs
				return(x)
			})

		raw <- map(raw, ~ .[score(.) >= threshold])

		raw_matrix <- select_samples %>%
			map(
				~as_tibble(., .name_repair = "unique") %>%
					mutate(position = paste(seqnames, start, end, strand, sep="_")) %>%
					select(position, score)
			) %>%
				bind_rows(.id = "sample") %>%
				spread(key = sample, value = score, fill = 0)
	}

	if (data_type == "tsr") {
		## Pull data from proper slot.
		if (samples == "all") samples <- names(experiment@experiment$TSRs)
		select_samples <- experiment@experiment$TSRs[samples]

		select_samples <- select_samples %>%
			map(function(x) {
				if ("nTAGs" %in% names(mcols(x))) score(x) <- x$nTAGs
				return(x)
			})

		select_samples <- select_samples %>%
			map(
				~ as_tibble(.x, .name_repair = "unique") %>%
					makeGRangesFromDataFrame(keep.extra.columns = TRUE)
			)

		raw <- map(select_samples, ~ .[score(.) >= threshold])

		## Merge overlapping TSRs to get consensus
		tsr_consensus <- experiment@experiment$TSRs %>%
			as("GRangesList") %>%
			unlist %>%
			reduce(ignore.strand=FALSE) %>%
			as_tibble(.name_repair = "unique") %>%
			mutate(names = paste(seqnames, start, end, strand, sep="_")) %>%
			makeGRangesFromDataFrame(keep.extra.columns = TRUE)

		raw_matrix <- map(
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
	}

	if (data_type == "features") {
		rnaseq_matrix <- experiment@experiment$features$rna_seq %>%
			rename(position = gene_id)

		fiveprime_matrix <- experiment@experiment$features$five_prime %>%
			rename(position = gene_id)

		raw_matrix <- left_join(rnaseq_matrix, fiveprime_matrix, by = "position")
	}

	## CPM normalize counts
	if (data_type %in% c("tss", "tsr")) {
		cpm_counts <- raw %>%
			map(function(x) {
				x$score <- x$score %>% cpm %>% as.numeric
				return(x)
			})
	}

	## Filter out positions that have less than n_samples number of samples with reads above threshold.
	filtered_matrix <- raw_matrix %>%
		mutate_if(is.numeric, ~ {.x >= threshold}) %>%
		mutate(rowsums = rowSums(.[, 2:ncol(.)])) %>%
		{which(.$rowsums >= n_samples)} %>%
		raw_matrix[., ]

	## TMM normalize counts.
	tmm_matrix <- filtered_matrix %>%
		column_to_rownames("position") %>%
		as.matrix %>%
		DGEList %>%
		calcNormFactors %>%
		cpm %>%
		as_tibble(.name_repair = "unique", rownames = "position")

	## Put counts into proper slots.
	if (data_type == "tss") {
		experiment@counts$TSSs <- list(
			"raw" = raw,
			"cpm" = cpm_counts,
			"raw_matrix" = filtered_matrix,
			"tmm_matrix" = tmm_matrix
		)
	} else if (data_type == "tsr") {
		experiment@counts$TSRs <- list(
			"raw" = raw,
			"cpm" = cpm_counts,
			"raw_matrix" = filtered_matrix,
			"tmm_matrix" = tmm_matrix
		)
	} else if (data_type == "features") {
		experiment@counts$features <- list(
			"raw_matrix" = filtered_matrix,
			"tmm_matrix" = tmm_matrix
		)
	}

	return(experiment)
}
