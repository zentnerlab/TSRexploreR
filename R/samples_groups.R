
#' Sample Sheet
#'
#' Add a sample sheet to a tsrexplorer object.
#'
#' @param experiment tsrexplorer object
#' @param sample_sheet Sample sheet as tab delimited file or data.frame.
#' Must have columns: 'tss_name' and 'tsr_name', and can have additional
#'   columns specifying condition, batch, etc.
#'
#' @rdname add_sample_sheet-function
#' @export

add_sample_sheet <- function(
  experiment,
  sample_sheet
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.string(sample_sheet) | is.data.frame(sample_sheet))

  ## Import sample sheet if it is a file.
  sample_sheet_type <- case_when(
    is.string(sample_sheet) ~ "character",
    is.data.frame(sample_sheet) ~ "data.frame"
  )

  sample_sheet <- switch(
    sample_sheet_type,
    "character"=fread(sample_sheet, sep="\t", header=TRUE),
    "data.frame"=sample_sheet
  )
  assert_that(has_name(sample_sheet, c("tss_name", "tsr_name")))

  ## Convert sample sheet to data.table
  if (!is.data.table(sample_sheet)) setDT(sample_sheet)
  
  ## Add sample sheet to tsrexplorer object.
  experiment@meta_data$sample_sheet <- sample_sheet

  return(experiment)
}

#' Feature Groups
#'
#' Add gene/transcript group info to a tsrexplorer object.
#'
#' @param experiment tsrexplorer object
#' @param group_sheet Group sheet as tab delimited file or data.frame
#' with column 'feature', and other named columns containing feature groups
#' @param feature_type either 'gene' or 'transcript'
#'
#' @rdname add_feature_group-function
#' @export

add_feature_group <- function(experiment, group_sheet, feature_type) {
  
  ## Import group sheet if it is a file.
  if (is(group_sheet, "character")) {
    group_sheet <- read.delim(
      group_sheet, sep="\t", header=TRUE,
      stringsAsFactors=FALSE
    )
  }

  ## Convert group sheet to data.table
  group_sheet <- as.data.table(group_sheet)

  ## Add group sheet to tsrexplorer object.
  experiment@meta_data$group_sheet <- group_sheet
  
  return(experiment)
}
