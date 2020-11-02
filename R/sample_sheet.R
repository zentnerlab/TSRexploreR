
#' Sample Sheet
#'
#' Add a sample sheet to a tsrexplorer object.
#'
#' @param experiment tsrexplorer object
#' @param sample_sheet Sample sheet as tab delimited file or data.frame.
#' Must have columns: 'sample_name' 'file_1' and 'file_2', and can have additional
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
  assert_that(
    is.data.frame(sample_sheet) || 
    (is.string(sample_sheet) && is.readable(sample_sheet))
  )

  ## Import sample sheet if it is a file.
  sample_sheet_type <- case_when(
    is.string(sample_sheet) ~ "file",
    is.data.frame(sample_sheet) ~ "data.frame"
  )

  sample_sheet <- switch(
    sample_sheet_type,
    "file"=fread(sample_sheet, sep="\t", header=TRUE),
    "data.frame"=as.data.table(sample_sheet)
  )

  assert_that(all(
    c("sample_name", "file_1", "file_2") %in%
    colnames(sample_sheet)
  ))

  ## Add sample sheet to tsrexplorer object.
  experiment@meta_data$sample_sheet <- sample_sheet

  return(experiment)
}
