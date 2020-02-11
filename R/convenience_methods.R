
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

extract_counts <- function(experiment, data_type, samples, cpm_norm = FALSE) {

	## Extract appropraite samples from TSSs or TSRs.
	if (data_type == "tss") {
		if (samples == "all") samples <- names(experiment@counts$TSSs$raw)
		selected_samples <- experiment@counts$TSSs$raw %>% extract(samples)
	} else if (data_type == "tsr") {
		if (samples == "all") samples <- names(experiment@counts$TSRs$raw)
		selected_samples <- experiment@counts$TSRs$raw %>% extract(samples)
	}

	## For each sample get the ranges and counts.
	counts <- selected_samples %>%
		map(function(x) {
			
			# Pull out the raw or cpm normalized counts.
			if (cpm_norm) {
				counts <- assay(x, "cpm")
			} else {
				counts <- assay(x, "raw")
			}
			counts <- counts %>%
				as_tibble(.name_repair = "unique") %>%
				rename(score = 1)

			# Add counts back to ranges.
			ranges <- rowRanges(x) %>%
				as_tibble(.name_repair = "unique") %>%
				bind_cols(counts)

			return(ranges)
		})

	## Return the ranges and counts as ouput of function.
	return(counts)
	
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

	## Get counts.
	if (data_type == "tss") {
		TSS_matrices
		if (samples == "all") samples <- names(experiment@counts$TSSs$normalized)
		selected_samples <- experiment@counts$TSSs$normalized %>%
			.[, .$sample %in% samples]
	} else if (data_type == "tsr") {
		if (samples == "all") samples <- names(experiment@counts$TSRs$normalized)
		selected_samples <- experiment@counts$TSRs$normalized %>%
			.[, .$sample %in% samples]
	}

	return(selected_samples)
} 
