#' tsrexplorer object
#'
#' @importFrom tibble tibble
#'
#' @slot experiment_info List of TSS GRanges returned by TSRchitect
#' @slot organism Organism name
#' @slot cores Number of CPU cores/threads available
#' @slot TMM TMM normalized TSS read counts

setClass(
	"tsr_object",
	representation(
		experiment_info = "list",
		organism = "character",
		cores = "numeric",
		TMM = "tbl"
	),
	prototype(
		experiment_info = list(),
		organism = NA_character_,
		cores = NA_real_,
		TMM = tibble()
	)
)

#' TSRexplorer constructor function.
#'
#' @param experiment_info List of TSS GRanges returned by TSRchitect
#' @param organism Name of organism associated with data
#' @param cores Number of CPU cores/threads available
#'
#' @return A tsrexplorer object
#'
#' @examples tsrexplorer(TSSs, "S. cerevisiae", 4)
#'
#' @importFrom GenomicRanges GRanges
#'
#' @export

tsrexplorer <- function(experiment_info, cores=1) {
		cores <- as.integer(cores)
		new("tsr_object", experiment_info=experiment_info, cores=cores)
}
