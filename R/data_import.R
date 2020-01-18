
#' Import TSSs
#'
#' Convenience function to import TSSs.
#'
#' @import tibble
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom dplyr rename mutate
#' @importFrom purrr map set_names discard
#' @importFrom magrittr %>%
#'
#' @param object Object with TSSs to import into tsrexplorer
#' @param tsrexplorer_obj tsrexplorer object to add the TSSs to
#' @param ... Additional arguments for classes
#'
#' @rdname tss_import-generic
#'
#' @export

setGeneric("tss_import", function(object, ...) {
	standardGeneric("tss_import")
})

#' Bedgraph/bigwig files sample sheet file
#'
#' @importFrom stringr str_detect regex
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "character"),
	function(tsrexplorer_obj, object) {
		## Prepare sample sheet.
		if (!file.exists(object)) {
			message(paste(object, "does not exist"))
			stop()
		} else {
			sample_sheet <- read.delim(object, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		}

		## Import data.
		imported_data <- sample_sheet %>%
			pmap(function(sample_name, pos, neg) {
				pos <- import(pos)
				neg <- import(neg)
				imported_data <- c(pos, neg)
			}) %>%
			set_names(pull(sample_sheet, sample_name))

		tss_experiment(tsrexplorer_obj) <- imported_data
		return(tsrexplorer_obj)
	}
)

#' Bedgraph/bigwig files data.frame
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "data.frame"),
	function(tsrexplorer_obj, object) {
		## Import data.
		imported_data <- sample_sheet %>%
			pmap(function(sample_name, pos, neg) {
				pos <- import(pos)
				neg <- import(neg)
				imported_data <- c(pos, neg)
			}) %>%
			set_names(pull(sample_sheet, sample_name))

		tss_experiment(tsrexplorer_obj) <- imported_data
		return(tsrexplorer_obj)
	}
)

#' TSRchitect object
#'
#' @importFrom TSRchitect tssObject
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "tssObject"),
	function(tsrexplorer_obj, object) {
		message("...Importing TSSs from TSRchitect object")
		imported_data <- object@tssCountData %>%
			rename(seqnames = seq, start = TSS, score = nTAGs) %>%
			mutate("end" = start) %>%
			makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
			set_names(object@sampleNames)

		tss_experiment(tsrexplorer_obj) <- imported_data
		return(tsrexplorer_obj)
	}
)

#' CAGEr object
#'
#' @importFrom CAGEr CAGEexp
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "CAGEexp"),
	function(tsrexplorer_obj, object) {
		message("Importing TSSs from CAGEexp object")
	}
)


