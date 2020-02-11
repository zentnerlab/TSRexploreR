#' TSRexplorer Class
#'
#' @slot experiment Named lists containing GRanges of TSSs and/or TSRs
#' @slot counts Named lists of TMM and CPM normalized TSSs and/or TSRs
#' @slot correlation Named lists of correlation values between TSS and/or TSR sets
#' @slot settings storage location for arguments used in various functions
#'
#' @rdname tsr_explorer-class
#' @export

setClass(
	"tsr_explorer",
	representation(
		experiment = "list",
		counts = "list",
		correlation = "list",
		settings = "list"
	),
	prototype(
		experiment = list(),
		counts = list(),
		correlation = list(),
		settings = list()
	)
)

#' TSRexplorer constructor function.
#'
#' @description
#' This function generates a new tsr_explorer object for
#' detailed analysis of transcription start sites (TSSs)
#' and TSS clusters, referred to as transcription
#' start regions (TSRs) or clustered transcription start sites (cTSSs).
#'
#' @import methods
#' @importFrom GenomicRanges GRanges
#' @importFrom tibble tibble
#'
#' @param TSSs Named list of TSS GRanges returned by TSRchitect
#' @param TSRs Named list of TSR GRanges returned by TSRchitect 
#'
#' @return A tsrexplorer object containing TSSs and/or TSRs
#'
#' @examples
#' TSSs <- system.file("extdata", "yeast_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#'
#' @rdname tsr_explorer-class
#' @export

tsr_explorer <- function(TSSs = NA, TSRs = NA) {

		tsr_obj <- new(
			"tsr_explorer",
			experiment = list("TSSs" = TSSs, "TSRs" = TSRs)
		)

		return(tsr_obj)
}
