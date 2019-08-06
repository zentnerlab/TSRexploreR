
setClass(
	"tsr_object",
	representation(
		experiment = "list",
		cores = "numeric",
		TMM = "tbl"
	),
	prototype(
		experiment = list(),
		cores = NA_real_,
		TMM = tibble()
	)
)

#' TSRexplorer constructor function.
#'
#' @param experiment List of TSS GRanges returned by TSRchitect
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

tsr_explorer <- function(experiment, cores=1) {
		cores <- as.integer(cores)
		new("tsr_object", experiment=experiment, cores=cores)
}
