#' TSRexplorer Class
#'
#' @slot experiment Named lists containing GRanges of TSSs and/or TSRs
#' @slot counts Named lists of TMM and CPM normalized TSSs and/or TSRs
#' @slot correlation Named lists of correlation values for TSS and/or TSR sets
#' @slot diff_features Differential features
#' @slot settings Storage location for arguments used in various functions
#' @slot meta_data Storage for meta_data (what metadata? qq)
#'
#' @rdname tsr_explorer-class
#'
#' @export

setClass(
  "tsr_explorer",
  representation(
    experiment="list",
    counts="list",
    correlation="list",
    diff_features="list",
    settings="list",
    meta_data="list"
  ),
  prototype(
    experiment=list(),
    counts=list(),
    correlation=list(),
    diff_features=list(),
    settings=list(),
    meta_data=list()
  )
)

#' TSRexplorer constructor function.
#'
#' @description
#' This function generates a new tsrexplorer object for
#' detailed analysis of transcription start sites (TSSs)
#' and TSS clusters, referred to here as transcription
#' start regions (TSRs).
#'
#' @import methods
#'
#' @param TSSs Named list of TSS GRanges
#' @param TSRs Named list of TSR GRanges
#'
#' @return A tsrexplorer object containing TSSs and/or TSRs
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#'
#' @rdname tsr_explorer-class
#'
#' @export

tsr_explorer <- function(TSSs=NA, TSRs=NA, sample_sheet=NULL) {

  ## Input Check.
  assert_that(
    is.na(TSSs) ||
   (is.list(TSSs) && has_attr(TSSs, "names"))
  )
  assert_that(
    is.na(TSRs) ||
    (is.list(TSRs) && has_attr(TSRs, "names"))
  )
  assert_that(
    is.null(sample_sheet) ||
    is.data.frame(sample_sheet) ||
    is.readable(sample_sheet)
  )

  ## Prepare sample sheet.
  sample_sheet_type <- case_when(
    is.null(sample_sheet) ~ "none",
    is.data.frame(sample_sheet) ~ "dataframe",
    is.character(sample_sheet) ~ "file"
  )

  if (sample_sheet_type != "none") {
    assert_that(
      c("sample_name", "file_1", "file_2") %in%
      colnames(sample_sheet)
    )
  }

  sample_sheet <- switch(
    sample_sheet_type,
    "none"=NA,
    "dataframe"=as.data.table(sample_sheet),
    "file"=fread(sample_sheet, sep="\t")
  )

  ## Create tsr explorer object.
  tsr_obj <- new(
    "tsr_explorer",
    experiment=list("TSSs"=TSSs, "TSRs"=TSRs),
    diff_features=list(
      "TSSs"=list(results=list()),
      "TSRs"=list(results=list())
    ),
    counts=list(
      "TSSs"=list(raw=list()),
      "TSRs"=list(raw=list())
    ),
    meta_data=list("sample_sheet"=sample_sheet)
  )

  return(tsr_obj)
}
