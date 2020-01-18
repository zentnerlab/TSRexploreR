
#' Import TSSs
#'
#' Convenience function to import TSSs.
#'
#' @import tibble
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom dplyr rename mutate
#' @importFrom TSRchitect tssObject
#' @importFrom purrr map set_names
#' @importFrom magrittr %>%
#' @importFrom CAGEr CAGEexp
#'
#' @param object Object with TSSs to import into tsrexplorer
#' @param ... Additional arguments for classes
#'
#' @rdname tss_import-generic
#'
#' @export

setGeneric("tss_import", function(object, ...) {
	standardGeneric("tss_import")
})

#' Directory of bedgraph/bigwig files
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "character"),
	function(object) {
		if (!is.vector(object) & dir.exists(object)) {
			data_files <- list.files(object, full.names = TRUE)

			imported_data <- data_files %>%
				map(~ import(.)) %>%
				set_names(basename(data_files))
		} else {
			imported_data <- object %>%
				map(~ import(.)) %>%
				set_names(basename(object))
		}

		tss_experiment(object) <- imported_data
		return(object)
	}
)

#' Data frame with single sample
#'
#' @param sample_name The name of the sample
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "data.frame"),
	function(object, sample_name = "sample") {
		## Check for required columns.
		required_columns <- c("seqnames", "start", "end", "strand", "score")

		if (!colnames(object) %>% {all(required_columns %in% .)}) {
			missing_columns <- discard(required_columns, ~ {.x %in% colnames(object)})
			message(paste("Data.frame missing columns:", paste0(missing_columns, collapse = " ")))
			stop()
		}

		converted_data <- object %>%
			makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
			list %>%
			set_names(sample_name)

		tss_experiment(object) <- converted_data
		return(object)
	}
)

#' TSRchitect object
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "tssObject"),
	function(object) {
		message("...Importing TSSs from TSRchitect object")
		imported_data <- object@tssCountData %>%
			rename(seqnames = seq, start = TSS, score = nTAGs) %>%
			mutate("end" = start) %>%
			makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
			set_names(object@sampleNames)

		tss_experiment(object) <- imported_data
		return(object)
	}
)

#' CAGEr object
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "CAGEexp"),
	function(object, data_type) {
		message("Importing TSSs from CAGEexp object")
	}
)
