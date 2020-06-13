
#' Add TSSs
#'
#' Convenience function to add TSSs to a tsrexplorer object.
#'
#' @include tsrexplorer.R
#'
#' @param tsrexplorer_object tsrexplorer object
#'
#' @rdname tss_experiment-generic
#'
#' @export

setGeneric("tss_experiment", function(tsrexplorer_object) standardGeneric("tss_experiment"))
setGeneric("tss_experiment<-", function(tsrexplorer_object, value) standardGeneric("tss_experiment<-"))

#' @rdname tss_experiment-generic

setMethod("tss_experiment", signature(tsrexplorer_object = "tsr_explorer"),
	function(tsrexplorer_object) {
		tsrexplorer_object@experiment$TSSs
	}
)

setMethod("tss_experiment<-", signature(tsrexplorer_object = "tsr_explorer"),
	function(tsrexplorer_object, value) {
		tsrexplorer_object@experiment$TSSs <- value
	}
)

#' Add TSRs
#'
#' Convenience function to add TSRs to a tsrexplorer object.
#'
#' @include tsrexplorer.R
#'
#' @param tsrexplorer_object tsrexplorer object
#'
#' @rdname tsr_experiment-generic
#'
#' @export

setGeneric("tsr_experiment", function(tsrexplorer_object) standardGeneric("tsr_experiment"))
setGeneric("tsr_experiment<-", function(tsrexplorer_object, value) standardGeneric("tsr_experiment<-"))

#' @rdname tsr_experiment-generic

setMethod("tsr_experiment", signature(tsrexplorer_object = "tsr_explorer"),
	function(tsrexplorer_object) {
		tsrexplorer_object@experiment$TSRs
	}
)

setMethod("tsr_experiment<-", signature(tsrexplorer_object = "tsr_explorer"),
        function(tsrexplorer_object, value) {
                tsrexplorer_object@experiment$TSRs <- value
        }
)

#' Extract Counts
#'
#' @description
#' Extract counts from a tsrexplorer object.
#'
#' @import tibble
#' @importFrom magrittr extract %>%
#' @importFrom SummarizedExperiment assay rowRanges
#' @importFrom purrr map
#' @importFrom dplyr bind_cols rename
#'
#' @param experiment tsrexplorer object
#' @param data_type Whether to extract TSS or TSR counts
#' @param samples Names of samples from which to extract counts
#' @param use_cpm Whether CPM values should be returned
#'
#' @rdname extract_counts-function
#' @export

extract_counts <- function(experiment, data_type, samples, use_cpm = FALSE) {

	## Extract appropriate TSS or TSR samples.
	if (data_type == "tss") {
		if (all(samples == "all")) samples <- names(experiment@counts$TSSs$raw)
		selected_samples <- experiment@counts$TSSs$raw[samples]
	} else if (data_type == "tsr") {
		if (all(samples == "all")) samples <- names(experiment@counts$TSRs$raw)
		selected_samples <- experiment@counts$TSRs$raw[samples]
	} else if (data_type == "tss_features") {
		if (all(samples == "all")) samples <- names(experiment@counts$TSS_features$raw)
		selected_samples <- experiment@counts$TSS_features$raw[samples]
	} else if (data_type == "tsr_features") {
		if (all(samples == "all")) samples <- names(experiment@counts$TSR_features$raw)
		selected_samples <- experiment@counts$TSR_features$raw[samples]
	}

	## Return counts as copies so you don't unintentionally overwrite the tsrexplorer object values.
	return_samples <- map(selected_samples, copy)

	## Return CPM score if requested.
	if (use_cpm) {
		walk(return_samples, function(x) {
			x[, score := cpm]
			x[, cpm := NULL]
			return(x)
		})
	}

	## Return the ranges and counts as output of function.
	return(return_samples)
	
}

#' Extract Count Matrices
#'
#' Extract normalized count matrices.
#'
#' @import tibble
#'
#' @param experiment tsrexplorer object
#' @param data_type Whether to extract 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param samples Samples to extract
#'
#' @rdname extract_matrix-function
#' @export

extract_matrix <- function(experiment, data_type, samples) {

	if (data_type == "tss") {
		selected_samples <- experiment@counts$TSSs$matrix
	} else if (data_type == "tsr") {
		selected_samples <- experiment@counts$TSRs$matrix
	} else if (data_type == "tss_features") {
		selected_samples <- experiment@counts$TSS_features$matrix
	} else if (data_type == "tsr_features") {
		selected_samples <- experiment@counts$TSR_features$matrix
	}

	if (samples == "all") samples <- colnames(selected_samples)
	selected_samples <- selected_samples[, samples]

	return(selected_samples)
}

#' Extract Differential Gene Sets
#'
#' Extract differential expression analysis results.
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The comparison sets to extract
#'
#' @rdname extract_de-function
#' @export

extract_de <- function(experiment, data_type, de_comparisons = "all") {
       	
	if (data_type == "tss") {
		if (de_comparisons == "all") de_comparisons <- names(experiment@diff_features$TSSs$results)
                de_samples <- experiment@diff_features$TSSs$results[de_comparisons]
		de <- map(de_samples, copy)
        } else if (data_type == "tsr") {
		if (de_comparisons == "all") de_comparisons <- names(experiment@diff_features$TSRs$results)
                de_samples <- experiment@diff_features$TSRs$results[de_comparisons]
		de <- map(de_samples, copy)
        } else if (data_type == "tss_features") {
		if (de_comparisons == "all") de_comparisons <- names(experiment@diff_features$TSS_features$results)
                de_samples <- experiment@diff_features$TSS_features$results[de_comparisons]
		de <- map(de_samples, copy)
        } else if (data_type == "tsr_features") {
		if (de_comparisons == "all") de_comparisons <- names(experiment@diff_features$TSR_features$results)
                de_samples <- experiment@diff_features$TSR_features$results[de_comparisons]
		de <- map(de_samples, copy)
        }

	return(de)
} 
