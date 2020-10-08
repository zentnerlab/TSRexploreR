
#' Import TSSs
#'
#' @description
#' Function to import TSSs from various sources.
#'
#' @param tsrexplorer tsrexplorer object
#' @param object Object with TSSs to import into tsrexplorer
#' @param ... Additional arguments for classes
#'
#' @details
#' TSRexploreR can import TSSs from various sources.
#' Currently bedgraphs, bigwigs, and TSRchitect tssObjects are supported.
#'
#' To import bedgraphs or bigwigs, a sample sheet must first be created.
#' The sample sheet can be a data.frame or tabular file.
#' It should have three columns: sample_name, file_1, and file_2.
#' sample_name specifies the sample name that will be added for that TSS in the
#'   tsr explorer object.
#' file_1 and file_2 should be the path to the two bedgraphs or bigwigs,
#'   which are usually separated by strand.
#'
#' To import TSSs directly from a TSRchitect tssObject,
#'   the TSRchitect workflow must first be run up to the
#'   'processTSS' step.
#'
#' @return tsr explorer object with imported TSSs
#'
#' @seealso
#' \code{\link{tsr_import}} to import TSRs.
#' \code{\link{tss_export}} to export TSSs.
#' \code{\link{tsr_export}} to export TSRs.
#'
#' @rdname tss_import-generic
#' @export

setGeneric(
  "tss_import",
  function(tsrexplorer_obj, object, ...) standardGeneric("tss_import"),
  signature = "object"
)

#' Bedgraph/bigwig files (sample sheet saved as file) qq sample sheet specifications?
#'
#' @param sep Delimiter for sample_sheet file
#' @param col_names Logical, whether the sample_sheet file has a header (TRUE) or not (FALSE)
#'
#' @rdname tss_import-generic
#'
#' @export

setMethod("tss_import", signature = signature(object = "character"),
  function(
    tsrexplorer_obj,
    object,
    sep = "\t",
    col_names = TRUE
  ) {

    ## Check inputs.
    if (!is(tsrexplorer_obj, "tsr_explorer")) stop("tsrexplorer_obj must be a tsr explorer object")

    if (!is(sep, "character") || length(sep) > 1) stop("sep must be a character")

    if (!is(col_names, "logical") || length(col_names) > 1) stop("col_names must be TRUE or FALSE")
  
    ## Check to see if sample sheet exists.
    if (!file.exists(object)) {
      message(str_c(object, "does not exist", sep = " "))
    }

    sample_sheet <- read.delim(object, sep = sep, header = col_names, stringsAsFactors = FALSE)

    ## Import data.
    imported_data <- sample_sheet %>%
      pmap(function(sample_name, file_1, file_2) {
        pos <- import(file_1)
        neg <- import(file_2)
        imported_data <- sort(c(pos, neg))
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
#'
#' @export

setMethod("tss_import", signature = signature(object = "data.frame"),
  function(
    tsrexplorer_obj,
    object
  ) {
    ## Check inputs.
    if (!is(tsrexplorer_obj, "tsr_explorer")) stop("tsrexplorer_obj must be a tsr explorer object")

    ## Import data.
    imported_data <- sample_sheet %>%
      pmap(function(sample_name, file_1, file_2) {
        pos <- import(file_1)
        neg <- import(file_2)
        imported_data <- sort(c(pos, neg))
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
  function(
    tsrexplorer_obj,
    object
  ) {
    ## Check inputs.
    if (!is(tsrexplorer_obj, "tsr_explorer")) stop("tsrexplorer_obj must be a tsr explorer object")

    ## Pull the TSSs out of the TSRchitect tssObject.
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
    
    ## Get TSR counts from CAGEexp object.
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
#' @description
#' Function to import TSRs from various sources.
#'
#' @param object Object with TSRs to import into tsrexplorer
#' @param tsrexplorer_obj tsrexplorer object to which to add TSRs
#' @param ... Additional arguments for classes
#'
#' @details
#' TSRexploreR can import TSRs from various sources.
#' Currently beds and TSRchitect tssObjects are supported.
#'
#' To import beds, a sample sheet must first be created.
#' The sample sheet can be a data.frame or tabular file.
#' It should have three columns: sample_name, file_1, and file_2.
#' sample_name specifies the sample name that will be added for that TSR in the
#'   tsr explorer object.
#' file_1 should be the path to the bed file,
#'   and file_2 can be left blank or filled with NA values.
#'
#' To import TSRs directly from a TSRchitect tssObject,
#'   the TSRchitect workflow must first be run up to the
#'   'determineTSR' step.
#'
#' @return tsr explorer object with added TSRs.
#'
#' @seealso
#' \code{\link{tss_import}} to import TSSs.
#' \code{\link{tss_export}} to export TSSs.
#' \code{\link{tsr_export}} to export TSRs.
#'
#' @rdname tsr_import-generic
#' @export

setGeneric(
  "tsr_import",
  function(tsrexplorer_obj, object, ...) standardGeneric("tsr_import"),
  signature = "object"
)

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
      pmap(function(sample_name, file_1, file_2) {
        imported_data <- import(file_1)
        return(imported_data)
      }) %>%
      set_names(pull(sample_sheet, "sample_name"))

    ## Add TSRs to tsrexplorer object.
    tsrexplorer_obj@experiment$TSRs <- imported_data
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
      pmap(function(sample_name, file_1, file_2) {
        imported_data <- import(file_1)
        return(imported_data)
      }) %>%
      set_names(pull(object, "sample_name"))

    ## Add TSRs to tsrexplorer object.
    tsrexplorer_obj@experiment$TSRs <- imported_data
    return(tsrexplorer_obj)
  }
)
