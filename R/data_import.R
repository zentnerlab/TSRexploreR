
#' Import TSSs and TSRs
#'
#' Convenience function to import TSSs and/or TSRs from various sources.
#'
#' @import tibble
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom dplyr rename mutate
#' @importFrom TSRchitect tssObject
#'
#' @param object Object to import into tsrexplorer
#' @param data_type Whether the data being imported is 'TSSs' or 'TSRs'
#'
#' @rdname tsrexplorer_import-generic
#'
#' @export

setGeneric("tsrexplorer_import", function(object, data_type) {
	standardGeneric("tsrexplorer_import")
})

#' Directory of bed/bedgraph/bigwig files
#'
#' @rdname tsrexplorer_import-generic

setMethod("tsrexplorer_import", signature(object = "character"),
	function(object, data_type) {
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
	}
)

#' Data frame with single sample
#'
#' @rdname tsrexplorer_import-generic

setMethod("tsrexplorer_import", signature(object = "data.frame"),
	function(object, data_type) {
		## Check for required columns.
		required_columns <- c("seqnames", "start", "end", "strand", "score")

		if (!colnames(object) %>% {all(required_columns %in% .)}) {
			missing_columns <- discard(required_columns, ~ {.x %in% colnames(object)})
			message(paste("Data.frame missing columns:", paste0(missing_columns, collapse = " ")))
			stop()
		}

		converted_data <- makeGRangesFromDataFrame(object, keep.extra.columns = TRUE)
	}
)

#' TSRchitect object
#'
#' @rdname tsrexplorer_import-generic

setMethod("tsrexplorer_import", signature(object = "tssObject"),
	function(object, data_type) {
		if (data_type == "TSSs") {
			message("importing TSSs from TSRchitect object")
			imported_data <- object@tssCountData %>%
				rename(seqnames = seq, start = TSS, score = nTAGs) %>%
				mutate("end" = start) %>%
				makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
				set_names(object@sampleNames)
		} else {
			message("importing TSRs from TSRchitect object")
		}
})
