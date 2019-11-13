
setClass(
	"tsr_object",
	representation(
		experiment = "list",
		samples = "vector",
		cores = "numeric",
		annotated = "list",
		counts = "list"
	),
	prototype(
		experiment = list(),
		samples = c(),
		cores = NA_real_,
		annotated = list(),
		counts = list()
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
#' @examples tsrexplorer(TSSs, 4)
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom tibble tibble
#'
#' @export

tsr_explorer <- function(TSSs, TSRs=NA, cores=1) {
		cores <- as.integer(cores)
		tsr_obj <- new(
			"tsr_object",
			samples=names(TSSs),
			cores=cores
		)

		tsr_obj@experiment$TSSs <- TSSs
		tsr_obj@experiment$TSRs <- TSRs

		return(tsr_obj)
}
