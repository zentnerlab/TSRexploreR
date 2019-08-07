
setClass(
	"tsr_object",
	representation(
		experiment = "list",
		samples = "vector",
		cores = "numeric",
		annotated = "list",
		normalized_counts = "list"
	),
	prototype(
		experiment = list(),
		samples = c(),
		cores = NA_real_,
		annotated = list(),
		normalized_counts = list()
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
		new("tsr_object", experiment=experiment, samples=names(experiment), cores=cores)
}
