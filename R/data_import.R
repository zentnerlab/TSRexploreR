
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
#' @importFrom stringr str_detect regex
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "character"),
	function(object) {
		## Check if string is a directory.
		if (!is.vector(object) & dir.exists(object)) {
			data_files <- list.files(object, full.names = TRUE)

			imported_data <- data_files %>%
				map(~ import(.)) %>%
				set_names(basename(data_files))
		## If string isn't a directory it's probably a file or vector of files.
		} else {
			data_files <- object
		}

		## Check if files have valid formats.
		if (!str_detect(data_files, regex("\\.(bedgraph|bigwig|bw|wig)$", ignore_case = TRUE)) %>% all) {
			error_files = data_files %>%
				discard(~ str_detect(., regex("\\.(bedgraph|bigwig|bw|wig)$", ignore_case = TRUE))) %>%
				paste0(collapse = " ")

			message(paste("Files not acceptable format:", error_files))
			stop()
		}

		## Import data.
		imported_data <- data_files %>%
			map(~ import(.)) %>%
			set_names(basename(data_files))

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
#' @importFrom TSRchitect tssObject
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
#' @importFrom CAGEr CAGEexp
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "CAGEexp"),
	function(object) {
		message("Importing TSSs from CAGEexp object")
	}
)


