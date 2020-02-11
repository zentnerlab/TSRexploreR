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
	data_type = c("tss", "tsr", "features"),
	samples = "all",
	threshold = 1,
	n_samples = 1
) {

	## Select proper slot of data.
	if (data_type == "tss") {
		count_matrix <- experiment@counts$TSSs$matrix
	} else if (data_type == "tsr") {
		count_matrix <- experiment@counts$TSRs$matrix
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

	row_ranges <- rowRanges(select_samples)
	col_data <- DataFrame(sample = colnames(select_samples))

	## Create filtered and TMM normalized RangedSummarizedExperiment.
	tmm_experiment <- SummarizedExperiment(
		assays = list("filtered" = filtered_counts, "tmm" = tmm_counts),
		rowRanges = row_ranges,
		colData = col_data
	)

	## Return the TMM normalized counts.
	experiment@correlation$tmm <- tmm_experiment
	return(experiment)
}
