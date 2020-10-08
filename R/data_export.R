
#' Export TSSs
#'
#' @description
#' Export TSSs to tables or bedgraphs
#'
#' @param experiment tsrexplorer object
#' @param samples Either 'all', or names of samples to export
#' @param file_type either 'bedgraph' or 'table'
#' @param out_dir Output directory for files
#' @param diff_tss If TRUE will output differential TSSs
#' @param sep Delimiter for table output
#'
#' @details
#' This function will save TSSs as bedgraphs, or a tab delimited file.
#'
#' 'file_type' controls what the TSSs will be output as.
#' 'bedgraph' will result in each sample being saved as two bedgraph files,
#'   one for each strand.
#' 'table' will output a file with the delimiter specified by the 'sep' argument.
#' The resulting table will have all columns added to the TSSs data
#'   in the tsr explorer object, such as various metrics.
#'
#' The directory to output the files can be set with 'out_dir'.
#'   A value of NA will save the files to the working directory by default.
#'
#' If 'diff_tss' is TRUE, only differential TSSs will be output.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tss_export(tsre_exp)
#'
#' @return Either bedgraphs split by strand, or a tabular file.
#'
#' @seealso
#' \code{\link{tsr_export}} to export TSRs.
#' \code{\link{tss_import}} to import TSSs.
#' \code{\link{tsr_import}} to import TSRs.
#'
#' @rdname tss_export-function
#' @export

tss_export <- function(
  experiment,
  samples = "all",
  file_type = "bedgraph",
  out_dir = NA,
  diff_tss = FALSE,
  sep = "\t"
) {

  ## Input checks.
  if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")

  if (!is(samples, "character")) stop("samples must be a character")

  if (!is(file_type, "character") || length(file_type) > 1) stop("file_type must be 'bedgraph' or 'table'")
  file_type <- str_to_lower(file_type)
  if (!file_type %in% c("bedgraph", "table")) stop("file_type must be 'bedgraph' or 'table'")

  if (!is.na(out_dir) && !is(out_dir, "character")) stop("out_dir must be NA or a character")

  if (!is(diff_tss, "logical")) stop("diff_tss must be TRUE or FALSE")

  ## Retrieve samples.
  if (!diff_tss) {
    if (all(samples == "all")) samples <- names(experiment@counts$TSSs$raw)
    export_samples <- experiment@counts$TSSs$raw[samples]
  } else {
    if (all(samples == "all")) samples <- names(experiment@diff_features$TSSs$results)
    export_samples <- experiment@diff_features$TSSs$results[samples]
  }

  ## Export files.
  if (file_type == "bedgraph") {
    iwalk(export_samples, function(x, y) {
      x <- sort(as_granges(x))
      
      pos_data <- x[strand(x) == "+"]
      pos_file <- file.path(
        ifelse(is.na(out_dir), getwd(), out_dir),
        str_c(y, "_pos.bedgraph")
      )
      export(pos_data, pos_file, "bedgraph")

      min_data <- x[strand(x) == "-"]
      min_file <- file.path(
        ifelse(is.na(out_dir), getwd(), out_dir),
        str_c(y, "_min.bedgraph")
      )
      export(min_data, min_file, "bedgraph")
    })
  } else if (data_type == "table") {
    iwalk(export_samples, function(x, y) {
      export_file <- file.path(
        ifelse(is.na(out_dir), getwd(), out_dir),
        str_c(y, "_TSSs.tsv")
      )
      fwrite(
        x, export_file, sep = sep, col.names = TRUE,
        row.names = FALSE, quote = FALSE
      )
    })
  }
}

#' Export TSRs
#'
#' @description
#' Export TSRs to table or BED
#'
#' @param experiment tsrexplorer object
#' @param samples Samples to export ('all' as well? qq)
#' @param file_type either 'bed' or 'table'
#' @param out_dir Output directory for files
#' @param diff_tsr Whether to pull out the differential TSRs (qq again, this is a bit unclear)
#' @param sep Delimiter for table output
#'
#' @details
#' This function will save TSRs as beds, or a tab delimited file.
#'
#' 'file_type' controls what the TSRs will be output as.
#' 'bed' will result in each sample being saved as a file.
#' 'table' will output a file with the delimiter specified by the 'sep' argument.
#' The resulting table will have all columns added to the TSRs data
#'   in the tsr explorer object, such as various metrics.
#'
#' The directory to output the files can be set with 'out_dir'.
#'   A value of NA will save the files to the working directory by default.
#'
#' If 'diff_tsr' is TRUE, only differential TSRs will be output.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#' tsr_export(tsre_exp)
#'
#' @return Either bed files, or a tabular file.
#'
#' @seealso
#' \code{\link{tss_export}} to export TSSs.
#' \code{\link{tss_import}} to import TSSs.
#' \code{\link{tsr_import}} to import TSRs.
#'
#' @rdname tsr_export-function
#' @export

tsr_export <- function(
  experiment,
  samples = "all",
  file_type = "bed",
  out_dir = NA,
  diff_tsr = FALSE,
  sep = "\t"
) {

  ## Check inputs.
  if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")

  if (!is(samples, "character")) stop("samples must be a character")

  if (!is(file_type, "character") || length(file_type) > 1) stop("file_type must be 'bed' or 'table'")
  file_type <- str_to_lower(file_type)
  if (!file_type %in% c("bed", "table")) stop("file_type must be 'bed' or 'table'")

  if (!is.na(out_dir) && !is(out_dir, "character")) stop("out_dir must be NA or a character")

  if (!is(diff_tsr, "logical")) stop("diff_tsr must be TRUE or FALSE")

  ## Retrieve samples.
  if (!diff_tsr) {
    if (all(samples == "all")) samples <- names(experiment@counts$TSRs$raw)
    export_samples <- experiment@counts$TSRs$raw[samples]
  } else {
    if (all(samples == "all")) samples <- names(experiment@diff_features$TSRs$results)
    export_samples <- experiment@diff_features$TSRs$results[samples]
  }

  ## Export files.
  if (file_type == "bed") {
    iwalk(export_samples, function(x, y) {
      x <- as_granges(x)

      bed_file <- file.path(
        ifelse(is.na(out_dir), getwd(), out_dir),
        str_c(y, "_TSRs.bed")
      )

      export(x, bed_file, "bed")
    })
  } else if (file_type == "table") {
    iwalk(export_samples, function(x, y) {
      export_file <- file.path(
        ifelse(is.na(out_dir), getwd(), out_dir),
        str_c(y, "_TSRs.tsv")
      )

      fwrite(
        x, export_file, sep = sep, col.names = TRUE,
        row.names = FALSE, quote = FALSE
      )
    })
  }

}
