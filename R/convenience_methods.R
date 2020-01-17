
#' Add TSSs or TSRs
#'
#' Convenience function to add TSSs or TSRs to tsrexplorer object
#'
#' @importClassesFrom tsrexplorer tsr_object
#'
#' @param fiveprime_data Named list of GRanges with TSS or TSR data
#' @param data_type Store either 'TSSs' or 'TSRs' into the tsrexplorer object
#'
#' @rdname fiveprime_data-generic
#'
#' @export

#setGeneric("fiveprime_data", function(object, data_type) standardGeneric("fiveprime_data"))
#setGeneric("fiveprime_data<-", function(object, data_type, fiveprime_data) standardGeneric("fiveprime_data<-"))

#' @rdname fiveprime_data-generic

#setMethod("fiveprime_data", signature(object = "tsr_object"),
#	function(object, data_type = c("TSSs", "TSRs")) {
#		if (data_type == "TSSs") {
#			tsrexplorer_object@experiment$TSSs
#		} else if (data_type == "TSRs") {
#			tsrexplorer_object@experiment$TSRs
#		}
#	}
#)

#setMethod("fiveprime_data<-", signature(object = "tsr_object"),
#	function(tsrexplorer_object, data_type = c("TSSs", "TSRs"), fiveprime_data) {
#		if (data_type == "TSSs") {
#			tsrexplorer_object@experiment$TSSs <- fiveprime_data
#		} else if (data_type == "TSRs") {
#			tsrexplorer_object@experiment$TSRs <- fiveprime_data
#		}
#	}
#)
