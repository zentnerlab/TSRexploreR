
#' Add TSSs
#'
#' Convenience function to add TSSs to tsrexplorer object
#'
#' @include tsrexplorer.R
#'
#' @param tsrexplorer_object tsr_explorer object
#'
#' @rdname tss_experiment-generic
#'
#' @export

setGeneric("tss_experiment", function(tsrexplorer_object) standardGeneric("tss_experiment"))
setGeneric("tss_experiment<-", function(tsrexplorer_object, fiveprime_data) standardGeneric("tss_experiment<-"))

#' @rdname tss_experiment-generic

setMethod("tss_experiment", signature(tsrexplorer_object = "tsr_explorer"),
	function(tsrexplorer_object) {
		tsrexplorer_object@experiment$TSSs
	}
)

setMethod("tss_experiment<-", signature(tsrexplorer_object = "tsr_explorer"),
	function(tsrexplorer_object, fiveprime_data) {
		tsrexplorer_object@experiment$TSSs <- fiveprime_data
	}
)

#' Add TSRs
#'
#' Convenience function to add TSRs to tsrexplorer object
#'
#' @include tsrexplorer.R
#'
#' @param tsrexplorer_object tsr_explorer object
#'
#' @rdname tsr_experiment-generic
#'
#' @export

setGeneric("tsr_experiment", function(tsrexplorer_object) standardGeneric("tsr_experiment"))
setGeneric("tsr_experiment<-", function(tsrexplorer_object, fiveprime_data) standardGeneric("tsr_experiment<-"))

#' @rdname tsr_experiment-generic

setMethod("tsr_experiment", signature(tsrexplorer_object = "tsr_explorer"),
	function(tsrexplorer_object) {
		tsrexplorer_object@experiment$TSRs
	}
)

setMethod("tsr_experiment<-", signature(tsrexplorer_object = "tsr_explorer"),
        function(tsrexplorer_object, fiveprime_data) {
                tsrexplorer_object@experiment$TSRs <- fiveprime_data
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
#' @param data_type whether to extract from the 'tss' or 'tsr' sets
#' @param samples names of samples to extract
#'
#' @rdname extract_counts-function
#' @export

extract_counts <- function(experiment, data_type, samples) {

	## Extract appropraite samples from TSSs or TSRs.
	if (data_type == "tss") {
		if (samples == "all") samples <- names(experiment@counts$TSSs$raw)
		selected_samples <- experiment@counts$TSSs$raw[samples]
	} else if (data_type == "tsr") {
		if (samples == "all") samples <- names(experiment@counts$TSRs$raw)
		selected_samples <- experiment@counts$TSRs$raw[samples]
	} else if (data_type == "tss_features") {
		if (samples == "all") samples <- names(experiment@counts$TSS_features$raw)
		selected_samples <- experiment@counts$TSS_features$raw[samples]
	} else if (data_type == "tsr_features") {
		if (samples == "all") samples <- names(experiment@counts$TSR_features$raw)
		selected_samples <- experiment@counts$TSR_features$raw[samples]
	}

	## Return the ranges and counts as ouput of function.
	return(selected_samples)
	
}

#' Extract Count Matrices
#'
#' Extract the normalized count matrices.
#'
#' @import tibble
#'
#' @param experiment tsrexplorer object
#' @param data_type Whether to extract 'tss', 'tsr', or 'feature' counts
#' @param samples Sampels to extract
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
#' Extract the differential expression results.
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The comparison sets to extract
#'
#' @rdname extract_de-function
#' @export

extract_de <- function(experiment, data_type, de_comparisons) {
       	
	if (data_type == "tss") {
                de_samples <- experiment@diff_features$TSSs
        } else if (data_type == "tsr") {
                de_samples <- experiment@diff_features$TSRs
        } else if (data_type == "tss_features") {
                de_samples <- experiment@diff_features$TSS_features
        } else if (data_type == "tsr_features") {
                de_samples <- experiment@diff_features$TSR_features
        }

        if (de_comparisons == "all") {
                de_samples <- discard(de_samples, names(de_samples) %in% c("model", "design"))
        } else {
                de_samples <- de_samples[de_comparisons]
        }

	return(de_samples)
} 
