#' CPM Normalize Counts
#'
#' CPM normalize the TSS, TSR, and/or feature counts.
#'
#' @param experiment tsrexplorer object
#' @param data_type 'tss', 'tsr', 'tss_features', 'tsr_features'
#'
#' @return tsrexplorer object with added CPM counts column.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- cpm_normalize(tsre_exp, data_type = "tss")
#'
#' @rdname cpm_normalize-function
#' @export

cpm_normalize <- function(experiment, data_type = c("tss", "tsr", "tss_features", "tsr_features")) {

	## Check inputs.
	if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")

	if (!is(data_type, "character")) stop("data_type must be a character")
	if (length(data_type) > 1) stop("data_type must be a character")
	data_type <- str_to_lower(data_type)
	if (!data_type %in% c("tss", "tsr", "tss_features", "tsr_features")) {
		stop("data_type must be 'tss', 'tsr', 'tss_features', or 'tsr_features'")
	}

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
			x <- x[, cpm := cpm(score)]
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
#' @description
#' Using edgeR to TMM normalize TSSs or TSRs.
#' This methods helps to reduce the influence that differing number of reads
#'   between samples can have on inter-sample comparisons.
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom SummarizedExperiment assay "assay<-"
#'
#' @param experiment tsrexplorer object
#' @param data_type Whether TSSs, TSRs, RNA-seq & five-prime feature counts should be normalized
#' @param threshold Filter out positions missing at least 'n_samples' number of samples with reads greater than or equal to threshold
#' @param n_samples Filter out positions missing at least n_samples number of samples with reads greater than or equal to 'threshold'
#' @param samples Vector with names of samples to include in the normalization
#'
#' @return tsr explorer object with tmm normalized count matrices
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- count_matrix(tsre_exp, data_type = "tss")
#' tsre_exp <- tmm_normalize(exp, data_type = "tss")
#'
#' @seealso \code{\link{count_matrix}} to prepare the matrices.
#'   \code{\link{plot_correlation}} for various correlation plots.
#'
#' @rdname tmm_normalize-function
#' @export

tmm_normalize <- function(
	experiment,
	data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	samples = "all",
	threshold = 1,
	n_samples = 1
) {

	## Check inputs.
	if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")

	if (!is(data_type, "character")) stop("data_type must be a character")
        if (length(data_type) > 1) stop("data_type must be a character")
        data_type <- str_to_lower(data_type)
        if (!data_type %in% c("tss", "tsr", "tss_features", "tsr_features")) {
                stop("data_type must be either 'tss', 'tsr', 'tss_features', or 'tsr_features'")
        }

	if (!is(samples, "character")) stop("samples must be a character vector")

	if (!is(threshold, "numeric") | !is(n_samples, "numeric")) {
		stop("threshold and n_samples must be positive integers")
	}
	if (threshold %% 1 != 0 | n_samples %% 1 != 0) {
		stop("threshold and n_samples must be positive integers")
	}
	if (threshold < 1 | n_samples < 1) {
		stop("threshold and n_samples must be greater than or equal to 1")
	}

	## Select proper slot of data.
	select_samples <- extract_matrix(experiment, data_type, samples)

	##  Filter the counts.
	sample_matrix <- as.data.table(assay(select_samples, "counts"))
	sample_matrix[, match := rowSums((.SD >= threshold)) >= n_samples]
	keep_ids <- which(sample_matrix[, match])
	
	select_samples <- select_samples[keep_ids, ]
	filtered_counts <- assay(select_samples, "counts")

	## TMM normalize filtered counts.
	tmm_counts <- filtered_counts %>%
		DGEList %>%
		calcNormFactors %>%
		cpm

	## Create filtered and TMM normalized RangedSummarizedExperiment.
	assay(select_samples, "tmm") <- tmm_counts

	## Return the TMM normalized counts.
	if (data_type == "tss") {
		experiment@counts$TSSs$matrix <- select_samples
	} else if (data_type == "tsr") {
		experiment@counts$TSRs$matrix <- select_samples
	} else if (data_type == "tss_features") {
		experiment@counts$TSS_features$matrix <- select_samples
	} else if (data_type == "tsr_features") {
		experiment@counts$TSR_features$matrix <- select_samples
	}

	return(experiment)
}
