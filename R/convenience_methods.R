
#' Add TSSs
#'
#' Convenience function to add TSSs to tsrexplorer object
#'
#' @include tsrexplorer.R
#'
#' @param tsrexplorer_object tsr_explorer object
#' @param fiveprime_data Named list of GRanges with TSSs
#'
#' @rdname tss_experiment-generic
#'
#' @export

setGeneric("tss_experiment", function(tsrexplorer_object, data_type) standardGeneric("tss_experiment"))
setGeneric("tss_experiment<-", function(tsrexplorer_object, data_type, fiveprime_data) standardGeneric("tss_experiment<-"))

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
