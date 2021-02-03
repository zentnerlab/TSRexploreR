#' Export TSSs
#'
#' @description
#' Export TSSs to bedGraph, bigWig, or table format.
#'
#' @inheritParams common_params
#' @param file_type Either 'bedgraph', 'table', or 'bigwig'.
#' @param out_dir Output directory for files.
#' @param diff_tss If TRUE, will output differential TSSs.
#' @param sep Delimiter for tabular output.
#'
#' @details
#' This function will save TSSs as bedGraphs, bigWigs, or a delimited table
#'
#' 'file_type' controls what the TSSs will be output as. 'bedgraph' will result 
#' in each sample being saved as two bedGraph files, one for each strand. The 
#' 'bigwig' argument provides a similar result. 'table' will output a file with 
#' the delimiter specified by the 'sep' argument. The resulting table will have 
#' all columns added to the TSSs data in the TSRexploreR object, such as annotation
#' information.
#'
#' The directory to output the files can be set with 'out_dir'. A value of NA will 
#' save the files to the working directory.
#'
#' If 'diff_tss' is TRUE, only differential TSSs will be output.
#'
#' @return Either bedGraphs or bigWigs split by strand, or a table.
#'
#' @seealso
#' \code{\link{tsr_export}} to export TSRs.
#' \code{\link{tss_import}} to import TSSs.
#' \code{\link{tsr_import}} to import TSRs.
#'
#' @examples
#' library("magrittr")
#' data(TSSs)
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss")
#' \donttest{tss_export(tsre)}
#'
#' @export

tss_export <- function(
  experiment,
  samples="all",
  file_type="bedgraph",
  out_dir=NA,
  diff_tss=FALSE,
  sep="\t",
  genome_assembly=NULL
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  file_type <- match.arg(
    str_to_lower(file_type),
    c("bedgraph", "table", "bigwig")
  )
  assert_that(is.na(out_dir) || is.character(out_dir))
  assert_that(is.flag(diff_tss))
  assert_that(
    is.null(genome_assembly) || is.character(genome_assembly) ||
    is(genome_assembly, "BSgenome")
  )

  ## Retrieve samples.
  if (!diff_tss) {
    if (all(samples == "all")) samples <- names(experiment@counts$TSSs$raw)
    export_samples <- experiment@counts$TSSs$raw[samples]
  } else {
    if (all(samples == "all")) samples <- names(experiment@diff_features$TSSs$results)
    export_samples <- experiment@diff_features$TSSs$results[samples]
  }

  ## If exporting as bigWig, retrieve the seq lengths from the genome object.
  genome_assembly <- .prepare_assembly(genome_assembly, experiment)

  assembly_type <- case_when(
    is(genome_assembly, "BSgenome") ~ "bsgenome",
    is(genome_assembly, "FaFile") ~ "fafile"
  )

  chrm_lengths <- switch(
    assembly_type,
    "fafile"=Rsamtools::seqinfo(genome_assembly),
    "bsgenome"=GenomeInfoDb::seqinfo(genome_assembly)
  )

  export_samples <- map(export_samples, function(x) {
    x <- sort(as_granges(x))
    chrm_lengths <- chrm_lengths[seqlevels(x)]
    seqlengths(x) <- seqlengths(chrm_lengths)
    return(x)
  })

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
  } else if (file_type == "table") {
    iwalk(export_samples, function(x, y) {
      export_file <- file.path(
        ifelse(is.na(out_dir), getwd(), out_dir),
        str_c(y, "_TSSs.tsv")
      )
      fwrite(
        x, export_file, sep=sep, col.names=TRUE,
        row.names=FALSE, quote=FALSE
      )
    })
  } else if (file_type == "bigwig") {
    iwalk(export_samples, function(x, y) {
      pos_data <- x[strand(x) == "+"]
      pos_file <- file.path(
        ifelse(is.na(out_dir), getwd(), out_dir),
        str_c(y, "_pos.bigwig") 
      )
      export(pos_data, pos_file, "bigwig")

      min_data <- x[strand(x) == "-"]
      min_file <- file.path(
        ifelse(is.na(out_dir), getwd(), out_dir),
        str_c(y, "_min.bigwig")
      )
      export(min_data, min_file, "bigwig")
    })
  }
}

#' Export TSRs
#'
#' @description
#' Export TSRs to table or BED
#'
#' @inheritParams common_params
#' @param file_type Either 'bed' or 'table'.
#' @param out_dir Output directory for files.
#' @param diff_tsr If TRUE, will output differential TSSs.
#' @param sep Delimiter for tabular output.
#'
#' @details
#' This function will save TSRs as BED files or a delimited table.
#'
#' 'file_type' controls what the TSRs will be output as. 'bed' will result in each 
#' sample being saved as a BED file. 'table' will output a file with the delimiter 
#' specified by the 'sep' argument. The resulting table will have all columns added 
#' to the TSR data in the TSRexplorer object, such as annotation information.
#'
#' The directory to output the files can be set with 'out_dir'. A value of NA will 
#' save the files to the working directory by default.
#'
#' If 'diff_tsr' is TRUE, only differential TSRs will be output.
#'
#' @return Either a BED file or a table.
#'
#' @seealso
#' \code{\link{tss_export}} to export TSSs.
#' \code{\link{tss_import}} to import TSSs.
#' \code{\link{tsr_import}} to import TSRs.
#'
#' @examples
#' library("magrittr")
#' data(TSSs)
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3)
#' \donttest{tsr_export(tsre)}
#'
#' @export

tsr_export <- function(
  experiment,
  samples="all",
  file_type="bed",
  out_dir=NA,
  diff_tsr=FALSE,
  sep="\t"
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  file_type <- match.arg(str_to_lower(file_type), c("bed", "table"))
  if (!is.na(out_dir) && !is(out_dir, "character")) stop("out_dir must be NA or a character")
  assert_that(is.flag(diff_tsr))

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
        x, export_file, sep=sep, col.names=TRUE,
        row.names=FALSE, quote=FALSE
      )
    })
  }

}
