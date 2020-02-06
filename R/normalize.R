#' TMM Normalize TSSs or TSRs
#'
#' Using edgeR to TMM normalize TSSs or TSRs
#'
#' @import tibble
#' @importFrom dplyr mutate mutate_at mutate_all vars select bind_rows group_by summarize mutate_if left_join rename
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom GenomicRanges makeGRangesFromDataFrame score "score<-" mcols "mcols<-"
#' @importFrom SummarizedExperiment SummarizedExperiment assay "assay<-"
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges findOverlapPairs
#' @importFrom tidyr spread complete
#' @importFrom purrr map imap
#' @importFrom magrittr %>% extract set_colnames
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

	# TSS data.
	if (data_type == "tss") {
		if (length(samples) == 1 & samples == "all") samples <- tss_experiment(experiment) %>% names
		select_samples <- tss_experiment(experiment) %>% extract(samples)

		raw_matrix <- select_samples %>%
			map(~ as_tibble(., .name_repair = "unique")) %>%
			bind_rows(.id = "sample") %>%
			spread(key = sample, value = score, fill = 0)
	}

	# TSR data.
	if (data_type == "tsr") {
		## Pull data from proper slot.
		if (length(samples) == 1 & samples == "all") samples <- tsr_experiment(experiment) %>% names
		select_samples <- tsr_experiment(experiment) %>% extract(samples)

		## Merge overlapping TSRs to get consensus
		tsr_consensus <- select_samples %>%
			purrr::reduce(c) %>%
			GenomicRanges::reduce(ignore.strand = FALSE)

		raw_matrix <- select_samples %>%
			map(
				~ findOverlapPairs(query = tsr_consensus, subject = .) %>%
					as_tibble(.name_repair = "unique") %>%
					select(first.seqnames, first.start, first.end, first.strand, second.X.score) %>%
					rename(
						seqnames = first.seqnames, start = first.start, end = first.end,
						strand = first.strand, score = second.X.score
					) %>%
					group_by(seqnames, start, end, strand) %>%
					summarize(score = sum(score)) %>%
					ungroup
			)  %>%
			bind_rows(.id = "sample") %>%
			spread(key = sample, value = score, fill = 0)
	}

	# Feature counts.
	if (data_type == "features") {
		rnaseq_matrix <- experiment@experiment$features$rna_seq
		fiveprime_matrix <- experiment@experiment$features$five_prime

		raw_matrix <- left_join(rnaseq_matrix, fiveprime_matrix, by = "gene_id")
	}

	## Store unfiltered and unfiltered-CPM normalized counts.
	if (data_type %in% c("tss", "tsr")) {
		raw_counts <- select_samples %>%
			imap(function(gr, sample_name) {
				count_data <- gr %>%
					score(.) %>%
					as.matrix %>%
					set_colnames(sample_name)
				cpm_data <- cpm(count_data)

				row_data <- gr
				score(row_data) <- NULL
				col_data <- DataFrame(sample = sample_name)

				raw_exp <- SummarizedExperiment(
					assays = list(raw = count_data, cpm = cpm_data),
					rowRanges = row_data,
					colData = col_data
				)
				return(raw_exp)
			})
	}

	## Store filtered, CPM, and TMM normalized count matrices.
	# Construct summarized experiment.
	if (data_type == "features") {
		count_matrix <- raw_matrix[,-1] %>% as.matrix
	} else {
		count_matrix <- raw_matrix[,-1:-5] %>% as.matrix
	}

	if (data_type == "features") {
		row_data <- DataFrame(raw_matrix[,1])
	} else {
		row_data <- makeGRangesFromDataFrame(raw_matrix)
	}

	col_data <- DataFrame("sample" = colnames(count_matrix), row.names = colnames(count_matrix))

	raw_exp <- SummarizedExperiment(
		assays = list(counts = count_matrix),
		rowRanges = row_data,
		colData = col_data
	)

	# Filter based on threshold and min samples.
	filtered_exp <- raw_exp %>%
		assay("counts") %>%
		as_tibble(.name_repair = "unique") %>%
		mutate_all(~ {.x >= threshold}) %>%
		mutate(rowsums = rowSums(.)) %>%
		{which(.$rowsums >= n_samples)} %>%
		raw_exp[., ]

	# CPM normalize counts.
	assay(filtered_exp, "cpm") <- filtered_exp %>%
		assay("counts") %>%
		cpm

	# TMM normalize counts.
	assay(filtered_exp, "tmm") <- filtered_exp %>%
		assay("counts") %>%
		DGEList %>%
		calcNormFactors %>%
		cpm

	## Put counts into proper slots.
	if (data_type == "tss") {
		experiment@counts$TSSs <- list(
			"raw" = raw_counts,
			"normalized" = filtered_exp
		)
	} else if (data_type == "tsr") {
		experiment@counts$TSRs <- list(
			"raw" = raw_counts,
			"normalized" = filtered_exp
		)
	} else if (data_type == "features") {
		experiment@counts$features <- list(
			"normalized" = filtered_exp
		)
	}

	return(experiment)
}
