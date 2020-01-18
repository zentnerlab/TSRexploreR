
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
