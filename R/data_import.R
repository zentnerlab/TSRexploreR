
#' Import TSSs
#'
#' Convenience function to import TSSs.
#'
#' @import tibble
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom dplyr rename mutate
#' @importFrom purrr pmap set_names
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

#' Bedgraph/bigwig files (sample sheet saved as file)
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

		## Add TSSs to tsrexplorer object.
		tss_experiment(tsrexplorer_obj) <- imported_data
		return(tsrexplorer_obj)
	}
)

#' Bedgraph/bigwig files (data.frame)
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

		## Add TSSs to tsrexplorer object.
		tss_experiment(tsrexplorer_obj) <- imported_data
		return(tsrexplorer_obj)
	}
)

#' TSRchitect object (tssObject)
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

		## Add TSSs to tsrexplorer object.
		tss_experiment(tsrexplorer_obj) <- imported_data
		return(tsrexplorer_obj)
	}
)

#' CAGEr object
#'
#' @importFrom CAGEr CAGEexp CTSStagCount
#'
#' @param data_type Either "tss", "tsr", or "consensus"
#'
#' @rdname tss_import-generic

setMethod("tss_import", signature(object = "CAGEexp"),
	function(tsrexplorer_obj, object, data_type) {
		message("Importing TSSs from CAGEexp object")

	if (data_type == "tss") {

		## Grab TSS counts from CAGEexp object.
		counts <- as.data.table(CTSStagCount(object))

		## Format counts for input to tsrexplorer.
		sample_names <- discard(colnames(counts), ~ . %in% c("chr", "pos", "strand"))
		counts <- melt(
			counts, variable.name = "sample", value.name = "score",
			measure.vars = sample_names, fill = 0
		)
		counts <- counts[score > 0]
		setnames(counts, old = c("chr", "pos"), new = c("seqnames", "start"))
		counts[, end := start]
		counts <- split(counts, counts$sample)

		## Change counts to GRanges objects.
		counts_gr <- counts %>%
			map(function(x) {
				x$sample <- NULL
				x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
				return(x)
			})

		## Add counts to tsrexplorer object.
		tss_experiment(tsrexplorer_obj) <- counts_gr

	} else if (data_type == "tsr") {
		
		## Grab TSR counts from CAGEexp object.
		counts <- tagClusters(object) %>%
			map(as.data.table)

		## Format counts for input to tsrexplorer.
		counts <- counts %>%
			map(function(x) {
				setnames(x, old = c("chr", "tpm"), new = c("seqnames", "score"))
				x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
				return(x)
			})

		## Add counts to tsrexplorer object.
		tsr_experiment(tsrexplorer_obj) <- counts
	
		return(tsrexplorer_obj)
	}
	}
)

#' Import TSRs
#'
#' Convenience function to import TSRs.
#'
#' @import tibble
#' @importFrom rtracklayer import
#' @importFrom purrr set_names pmap
#' @importFrom magrittr %>%
#'
#' @param object Object with TSRs to import into tsrexplorer
#' @param tsrexplorer_obj tsrexplorer object to add the TSRs to
#' @param ... Additional arguments for classes
#'
#' @rdname tsr_import-generic
#'
#' @export

setGeneric("tsr_import", function(object, ...) {
        standardGeneric("tsr_import")
})

#' Bed files (sample sheet saved as file) 
#'
#' @rdname tsr_import-generic

setMethod("tsr_import", signature(object = "character"),
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
			pmap(function(sample_name, bed) {imported_data <- import(bed)}) %>%
			set_names(pull(sample_sheet, sample_name))

		## Add TSRs to tsrexplorer object.
		tsrexplorer_obj <- tsr_experiment(imported_data)
		return(tsrexplorer_obj)
	}
)

#' Bed files (data.frame)
#'
#' @rdname tsr_import-generic

setMethod("tsr_import", signature(object = "data.frame"),
	function(tsrexplorer_obj, object) {
		## Import data.
		imported_data <- object %>%
			pmap(function(sample_name, bed) {imported_data <- import(bed)}) %>%
			set_names(pull(object, sample_name))

		## Add TSRs to tsrexplorer object.
		tsrexplorer_obj <- tsr_experiment(imported_data)
		return(tsrexplorer_obj)
	}
)
