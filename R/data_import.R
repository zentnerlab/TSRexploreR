#' Import TSSs
#'
#' @description
#' Function to import TSSs from various sources.
#'
#' @param experiment TSRexploreR object
#' @param sample_sheet Sample sheet
#' @param file_type Either 'auto', 'bedgraph', 'bigwig', or 'table'
#' @param delim If the input is a table, use this delimiter
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
#' @export

tss_import <- function(
  experiment,
  sample_sheet=NULL,
  file_type="auto",
  delim="\t"
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is.null(sample_sheet) ||
    (is.character(sample_sheet) | is.data.frame(sample_sheet))
  )
  file_type <- match.arg(
    str_to_lower(file_type),
    c("auto", "bigwig", "bedgraph", "table")
  )
  assert_that(is.string(delim))

  ## Convert sample sheet to data.table.
  sheet_type <- case_when(
    is.null(sample_sheet) & !is.null(experiment@meta_data$sample_sheet) ~ "internal",
    is.null(sample_sheet) & is.null(experiment@meta_data$sample_sheet) ~ "none",
    is.character(sample_sheet) ~ "file",
    is.data.frame(sample_sheet) ~ "table"
  )

  assert_that(sheet_type != "none")

  sample_sheet <- switch(
    sheet_type,
    "internal"=experiment@meta_data$sample_sheet,
    "file"=fread(sample_sheet, sep="\t", header=TRUE),
    "table"=as.data.table(sample_sheet)
  )

  ## Try to figure out file type if not specified.
  if (file_type == "auto") {
    file_ext <- sample_sheet[, .(file_1, file_2)] %>%
      unlist %>%
      discard(is.na) %>%
      str_extract("(?<=\\.)[[:alpha:]]+$") %>%
      str_to_lower %>%
      unique

    assert_that(
      is.string(file_ext),
      msg="All files must have the same file extension."
    )

    file_type <- case_when(
      file_ext %in% c("bw", "bigwig") ~ "bigwig",
      file_ext %in% c("csv", "tsv", "txt") ~ "table",
      file_ext == "bedgraph" ~ "bedgraph",
      TRUE ~ "unknown"
    )

    assert_that(
      file_ext != "unknown",
      msg=str_c(file_ext, "is an unknown format", sep=" ")
    )
  }

  ## Import TSSs.
  TSSs <- switch(
    file_type,
    "bedgraph"=.import_bedgraphs(sample_sheet),
    "bigwig"=.import_bigwigs(sample_sheet),
    "table"=.import_tables(sample_sheet, delim)
  )
  
  ## Add TSSs to TSRexploreR object.
  experiment@experiment$TSSs <- TSSs

  ## Add the sample sheet if provided to the function.
  if (!is.null(sample_sheet)) {
    experiment@meta_data$sample_sheet <- sample_sheet
  }

  return(experiment)
}

#' Import tables
#'
#' @param sample_sheet Sample sheet

.import_tables <- function(sample_sheet, delim) {
  samples <- sample_sheet[, .(sample_name, file_1)] %>%
    split(by="sample_name", keep.by=FALSE) %>%
    map(function(x) {
      x <- x %>%
        as.character %>%
        fread(sep=delim) %>%
        as_granges %>%
        sort
      return(x)
    })

  return(samples)
}

#' Import TSS Bedgraphs
#'
#' @param sample_sheet Sample sheet

.import_bedgraphs <- function(sample_sheet) {
  samples <- sample_sheet[, .(sample_name, file_1, file_2)] %>%
    split(by="sample_name", keep.by=FALSE) %>%
    map(function(x) {
      # Import positive strand.
      pos=x[, file_1]
      pos <- import(pos, "bedgraph")
      strand(pos) <- "+"

      # Import minus strand.
      neg=x[, file_2]
      neg <- import(neg, "bedgraph")
      strand(neg) <- "-"

      # Combine positive and negative strand.
      combined <- bind_ranges(pos, neg)

      # Sort the ranges.
      combined <- sort(combined)

      return(combined)
    })

  return(samples)
}

#' Import TSS bigwigs.
#'
#' @param sample_sheet Sample sheet

.import_bigwigs <- function(sample_sheet) {
  samples <- sample_sheet[, .(sample_name, file_1, file_2)] %>%
    split(by="sample_name", keep.by=FALSE) %>%
    map(function(x) {
      # Import positive strand.
      pos=x[, file_1]
      pos <- import(pos, "bigwig")
      strand(pos) <- "+"

      # Import minus strand.
      neg=x[, file_2]
      neg <- import(neg, "bigwig")
      strand(neg) <- "-"

      # Combine positive and negative strand.
      combined <- bind_ranges(pos, neg)

      # Sort the ranges.
      combined <- sort(combined)

      return(combined)
    })

  return(samples)
}


#' Import TSRs
#'
#' @description
#' Function to import TSRs from various sources.
#'
#' @param experiment TSRexploreR object
#' @param sample_sheet Sample sheet
#' @param file_type Either 'auto', 'table', or 'bed'
#' @param delim If input is a table specify the delimiter
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
#' @export

tsr_import <- function(
  experiment,
  sample_sheet=NULL,
  file_type="auto",
  delim="\t"
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is.null(sample_sheet) ||
    (is.character(sample_sheet) | is.data.frame(sample_sheet))
  )
  file_type <- match.arg(
    str_to_lower(file_type),
    c("auto", "table", "bed")
  )
  assert_that(is.string(delim))

  ## Convert sample sheet to data.table.
  sheet_type <- case_when(
    is.null(sample_sheet) & !is.null(experiment@meta_data$sample_sheet) ~ "internal",
    is.null(sample_sheet) & is.null(experiment@meta_data$sample_sheet) ~ "none",
    is.character(sample_sheet) ~ "file",
    is.data.frame(sample_sheet) ~ "table"
  )

  assert_that(sheet_type != "none")

  sample_sheet <- switch(
    sheet_type,
    "internal"=experiment@meta_data$sample_sheet,
    "file"=fread(sample_sheet, sep="\t", header=TRUE),
    "table"=as.data.table(sample_sheet)
  )

  ## Try to figure out file type if not specified.
  if (file_type == "auto") {
    file_ext <- sample_sheet[, .(file_1, file_2)] %>%
      unlist %>%
      discard(is.na) %>%
      str_extract("(?<=\\.)[[:alpha:]]+$") %>%
      str_to_lower %>%
      unique

    assert_that(
      is.string(file_ext),
      msg="All files must have the same file extension."
    )

    file_type <- case_when(
      file_ext %in% c("csv", "tsv", "txt") ~ "table",
      file_ext == "bed" ~ "bed",
      TRUE ~ "unknown"
    )

    assert_that(
      file_ext != "unknown",
      msg=str_c(file_ext, "is an unknown format", sep=" ")
    )
  }

  ## Import TSRs.
  TSRs <- switch(
    file_type,
    "bed"=.import_beds(sample_sheet),
    "table"=.import_tables(sample_sheet, delim)
  )

  ## Add TSRs to TSRexploreR object.
  experiment@experiment$TSRs <- TSRs

  ## Store in the sample sheet in the object if provided.
  if (!is.null(sample_sheet)) {
    experiment@meta_data$sample_sheet <- sample_sheet
  }

  return(experiment)

}

#' Import BEDs.
#'
#' @param sample_sheet Sample sheet

.import_beds <- function(sample_sheet) {

  samples <- sample_sheet[, .(sample_name, file_1)] %>%
    split(by="sample_name", keep.by=FALSE) %>%
    map(function(x) {
      x <- x %>%
        as.character %>%
        import("bed") %>%
        sort
      return(x)
    })

  return(samples)
}
