#' TSRexplorer Class
#'
#' @slot experiment Named lists containing GRanges of TSSs and/or TSRs
#' @slot cores Integer value for number of cores available
#' @slot annotated Named lists of annotated TSSs and/or TSRs
#' @slot counts Named lists of TMM and CPM normalized TSSs and/or TSRs
#' @slot correlation Named lists of correlation values between TSS and/or TSR sets
#'
#' @rdname tsr_explorer-class
#'
#' @export

setClass(
	"tsr_explorer",
	representation(
		experiment = "list",
		cores = "numeric",
		annotated = "list",
		counts = "list",
		correlation = "list"
	),
	prototype(
		experiment = list(),
		cores = NA_real_,
		annotated = list(),
		counts = list(),
		correlation = list()
	)
)

#' TSRexplorer constructor function.
#'
#' @param TSSs Named list of TSS GRanges returned by TSRchitect
#' @param TSRs Named list of TSR GRanges returned by TSRchitect 
#' @param cores Number of CPU cores/threads available
#'
#' @return A tsrexplorer object
#'
#' @import methods
#' @importFrom GenomicRanges GRanges
#' @importFrom tibble tibble
#'
#' @rdname tsr_explorer-class
#'
#' @export

tsr_explorer <- function(TSSs = NA, TSRs = NA, cores = 1) {

		tsr_obj <- new(
			"tsr_object",
			experiment = list("TSSs" = TSSs, "TSRs" = TSRs),
			cores = cores
		)

		return(tsr_obj)
}
