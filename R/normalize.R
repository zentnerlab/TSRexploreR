#' CPM Normalize Counts
#'
#' CPM normalize the TSS, TSR, and/or feature counts.
#'
#' @param experiment tsrexplorer object
#' @param data_type 'tss', 'tsr', 'tss_features', 'tsr_features'
#'
#' @rdname cpm_normalize-function
#' @export

cpm_normalize <- function(experiment, data_type = c("tss", "tsr", "tss_features", "tsr_features")) {
	
	## Grab appropriate samples.
	if (data_type == "tss") {
		select_samples <- experiment@counts$TSSs$raw
	} else if (data_type == "tsr") {
		select_samples <- experiment@counts$TSRs$raw
	} else if (data_type == "tss_features") {
		select_samples <- experiment@counts$TSS_features$raw
	} else if (data_type == "tsr_features") {
		select_samples <- experiment@counts$TSR_features$raw
	}

	## CPM normalize counts.
	cpm_counts <- select_samples %>%
		map(function(x) {
			cpm_matrix <- cpm(assay(x, "raw"))
			assay(x, "cpm") <- cpm_matrix
			return(x)
		})

	## Add data back to object.
	if (data_type == "tss") {
		experiment@counts$TSSs$raw <- cpm_counts
	} else if (data_type == "tsr") {
		experiment@counts$TSRs$raw <- cpm_counts
	} else if (data_type == "tss_features") {
		experiment@counts$TSS_features$raw <- cpm_counts
	} else if (data_type == "tsr_features") {
		experiment@counts$TSR_features$raw <- cpm_counts
	}

	return(experiment)
}

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
#' @rdname tmm_normalize-function
#'
#' @export

tmm_normalize <- function(
	experiment,
	data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	samples = "all",
	threshold = 1,
	n_samples = 1
) {

	## Select proper slot of data.
	if (data_type == "tss") {
		count_matrix <- experiment@counts$TSSs$matrix
	} else if (data_type == "tsr") {
		count_matrix <- experiment@counts$TSRs$matrix
	} else if (data_type == "tss_features") {
		count_matrix <- experiment@counts$TSS_features$matrix
	} else if (data_type == "tsr_features") {
		count_matrix <- experiment@counts$TSR_features$matrix
	}

	if (samples == "all") samples <- colnames(count_matrix)
	select_samples <- count_matrix[, count_matrix$sample %in% samples]

	##  Filter the counts.
	sample_matrix <- as.data.table(assay(select_samples, "counts"))
	sample_matrix[, match := rowSums((.SD >= threshold)) >= n_samples]
	keep_ids <- which(sample_matrix[, match])
	
	select_samples[keep_ids,]
	filtered_counts <- assay(select_samples, "counts")

	## TMM normalize filtered counts.
	tmm_counts <- select_samples %>%
		assay("counts") %>%
		DGEList %>%
		calcNormFactors %>%
		cpm

	if (data_type %in% c("tss_features", "tsr_features")) {
		row_data <- DataFrame("feature" = rownames(select_samples))
	} else {
		row_ranges <- rowRanges(select_samples)
	}
	col_data <- DataFrame(sample = colnames(select_samples))

	## Create filtered and TMM normalized RangedSummarizedExperiment.
	if (data_type %in% c("tss_features", "tsr_features")) {
		tmm_experiment <- SummarizedExperiment(
			assays = list("filtered" = filtered_counts, "tmm" = tmm_counts),
			rowData = row_data,
			colData = col_data
		)
	} else {
		tmm_experiment <- SummarizedExperiment(
			assays = list("filtered" = filtered_counts, "tmm" = tmm_counts),
			rowRanges = row_ranges,
			colData = col_data
		)
	}

	## Return the TMM normalized counts.
	if (data_type == "tss") {
		experiment@correlation$TSSs$tmm <- tmm_experiment
	} else if (data_type == "tsr") {
		experiment@correlation$TSRs$tmm <- tmm_experiment
	} else if (data_type == "tss_features") {
		experiment@correlation$TSS_features$tmm <- tmm_experiment
	} else if (data_type == "tsr_features") {
		experiment@correlation$TSR_features$tmm <- tmm_experiment
	}

	return(experiment)
}
