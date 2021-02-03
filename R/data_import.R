#' Import TSSs
#'
#' @description
#' Function to import TSSs from various sources.
#'
#' @inheritParams common_params
#' @param file_type Either 'auto', 'bedgraph', 'bigwig', 'table', or 'ctss'.
#' @param delim Delimiter for tabular input.
#'
#' @details
#' TSRexploreR can import TSSs from various file formats. Currently, these are bedGraphs, 
#' bigWigs, CTSS files, delimited tables, and TSRchitect tssObjects.
#'
#' To import TSSs from files, a sample sheet must first be created. The sample 
#' sheet can be a data.frame or tabular file. It should have three columns: 
#' sample_name, file_1, and file_2. sample_name specifies the sample name that 
#' will be added for that TSS in the TSRexploreR object. For bedGraphs and bigWigs,
#' file_1 and file_2 should be the paths to the plus and minus strand TSS files. For
#' CTSS or tabular import, file_1 should be the path to the file, and file_2can be 
#' left blank or filled with NA values.
#'
#' To import TSSs directly from a TSRchitect tssObject, the TSRchitect workflow 
#' must first be run up to the 'processTSS' step.
#'
#' @return TSRexploreR object with imported TSSs.
#'
#' @seealso
#' \code{\link{tsr_import}} to import TSRs.
#' \code{\link{tss_export}} to export TSSs.
#' \code{\link{tsr_export}} to export TSRs.
#'
#' @examples
#' \donttest{
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#'
#' # Export bedgraphs as example data.
#' tsre <- TSSs[1] %>%
#'   tsr_explorer %>%
#'   format_count(data_type="tss")
#' tss_export(tsre)
#'
#' # Import the previously created bedgraphs.
#' samples <- data.frame(
#'   sample_name="S288C_D_1",
#'   file_1="S288C_D_1_pos.bedgraph",
#'   file_2="S288C_D_1_min.bedgraph"
#' )
#' tss <- tsr_explorer(sample_sheet=samples)
#' tss_import(tss)
#' }
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
    c("auto", "bigwig", "bedgraph", "table", "ctss")
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
      file_ext == "ctss" ~ "ctss",
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
    "table"=.import_tables(sample_sheet, delim),
    "ctss"=.import_ctss(sample_sheet)
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
#' @inheritParams common_params
#' @param delim Delimiter for table.

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
#' @inheritParams common_params

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
#' @inheritParams common_params

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

#' Import CTSSs
#'
#' @inheritParams common_params

.import_ctss <- function(sample_sheet) {
  samples <- sample_sheet[, .(sample_name, file_1)] %>%
    split(by="sample_name", keep.by=FALSE) %>%
    map(function(x) {
      x <- x %>%
        as.character %>%
        fread(sep="\t", header=FALSE)
      setnames(
        x, old=seq_len(4),
        new=c("seqnames", "start", "strand", "score")
      )
      x[, end:=start]
      x <- sort(as_granges(x))

      return(x)
    })

  return(samples)
}

#' Import TSRs
#'
#' @description
#' Function to import TSRs from various sources.
#'
#' @inheritParams common_params
#' @param file_type Either 'auto', 'table', or 'bed'.
#' @param delim Delimiter for tabular input.
#'
#' @details
#' TSRexploreR can import TSRs from various sources. Currently, these are BED files,
#' delimited tables, and TSRchitect tssObjects.
#'
#' To import TSRs from files, a sample sheet must first be created. The sample sheet can 
#' be a data.frame or tabular file. It should have three columns: sample_name, file_1, 
#' and file_2. sample_name specifies the sample name that will be added for that 
#' TSR in the TSRexploreR object. file_1 should be the path to the file, and 
#' file_2 can be left blank or filled with NA values.
#'
#' To import TSRs directly from a TSRchitect tssObject, the TSRchitect workflow must 
#' first be run up to the 'determineTSR' step.
#'
#' @return TSRexploreR object with added TSRs.
#'
#' @seealso
#' \code{\link{tss_import}} to import TSSs.
#' \code{\link{tss_export}} to export TSSs.
#' \code{\link{tsr_export}} to export TSRs.
#'
#' @examples
#' \donttest{
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#'
#' # Export bed as example data.
#' tsre <- TSSs[1] %>%
#'   tsr_explorer %>%
#'   format_count(data_type="tss") %>%
#'   tss_clustering(threshold=3)
#' tsr_export(tsre)
#'
#' # Import the previously created bed.
#' samples <- data.frame(
#'   sample_name="S288C_D_1",
#'   file_1="S288C_D_1.bed",
#'   file_2=NA
#' )
#' tsr <- tsr_explorer(sample_sheet=samples)
#' tsr_import(tsr)
#' }
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
#' @inheritParams common_params

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
