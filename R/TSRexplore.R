#' TSRexploreR Class
#'
#' @slot experiment Named lists containing GRanges of TSSs and/or TSRs
#' @slot counts Named lists of TMM and CPM normalized TSSs and/or TSRs
#' @slot correlation Named lists of correlation values for TSS and/or TSR sets
#' @slot diff_features Differential features
#' @slot shifting TSS shifting data
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
    shifting="list",
    settings="list",
    meta_data="list"
  ),
  prototype(
    experiment=list(),
    counts=list(),
    correlation=list(),
    diff_features=list(),
    shifting=list(),
    settings=list(),
    meta_data=list()
  )
)

#' TSRexploreR constructor function.
#'
#' @description
#' This function generates a new TSRexploreR object for
#' detailed analysis of transcription start sites (TSSs)
#' and TSS clusters, referred to here as transcription
#' start regions (TSRs).
#'
#' @import methods
#'
#' @param TSSs Named list of TSS GRanges
#' @param TSRs Named list of TSR GRanges
#' @param sample_sheet Sample sheet
#' @param genome_annotation Genome annotation
#' @param genome_assembly Genome assembly
#'
#' @return A TSRexploreR object containing TSSs and/or TSRs
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#'
#' @rdname tsr_explorer-class
#'
#' @export

tsr_explorer <- function(
  TSSs=NA,
  TSRs=NA,
  sample_sheet=NULL,
  genome_annotation=NULL,
  genome_assembly=NULL
) {

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
    is.character(sample_sheet)
  )

  ## Prepare sample sheet.
  sample_sheet_type <- case_when(
    is.null(sample_sheet) ~ "none",
    is.data.frame(sample_sheet) ~ "dataframe",
    is.character(sample_sheet) ~ "file"
  )

  sample_sheet <- switch(
    sample_sheet_type,
    "none"=NA,
    "dataframe"=as.data.table(sample_sheet),
    "file"=fread(sample_sheet, sep="\t")
  )

  if (sample_sheet_type != "none") {
    assert_that(
      all(c("sample_name", "file_1", "file_2") %in%
      colnames(sample_sheet))
    )
  }

  ## Prepare genome assembly.
  genome_assembly <- .prepare_assembly(genome_assembly)

  ## Prepare genome annotation.
  genome_annotation <- .prepare_annotation(genome_annotation)

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
    shifting=list(results=list()),
    meta_data=list(
      "sample_sheet"=sample_sheet,
      "genome_annotation"=genome_annotation,
      "genome_assembly"=genome_assembly
    )
  )

  return(tsr_obj)
}

#' Prepare Genome Annotation
#'
#' @param annotation Genome annotation
#' @param experiment tsr explorer object

.prepare_annotation <- function(annotation, experiment=NULL) {

  ## Get annotation type.
  anno_type <- case_when(
    is.null(annotation) & !is.null(experiment) ~ "internal",
    is.null(annotation) & is.null(experiment) ~ "none",
    is.character(annotation) ~ "file",
    is(annotation, "TxDb") ~ "txdb"
  )

  if(anno_type == "internal") {
    assert_that(!is.null(experiment@meta_data$genome_annotation))
  }

  ## Prepare annotation based on type.
  annotation <- switch(
    anno_type,
    "internal"=experiment@meta_data$genome_annotation,
    "none"=NULL,
    "file"=makeTxDbFromGFF(annotation),
    "txdb"=annotation
  )

  return(annotation)

}

#' Prepare Genome Assembly
#'
#' @param assembly Genome assembly
#' @param experiment tsr explorer object

.prepare_assembly <- function(assembly, experiment=NULL) {

  ## Set the assembly type based on input.
  assembly_type <- case_when(
    is.null(assembly) & !is.null(experiment) ~ "internal",
    is.null(assembly) & is.null(experiment) ~ "none",
    is.character(assembly) ~ "file",
    is(assembly, "BSgenome") ~ "bsgenome"
  )

  ## Make sure an assembly is in the object if none provided.
  if (assembly_type == "internal") {
    assert_that(!is.null(experiment@meta_data$genome_assembly))
  }

  ## If assembly is a fasta file make sure it's indexed.
  if (
    assembly_type == "file" &&
    !file.exists(str_c(assembly, ".fai"))
  ) {
    indexFa(assembly)
  }

  ## Return the appropriate assembly.
  assembly <-switch(
    assembly_type,
    "none"=NULL,
    "internal"=experiment@meta_data$genome_assembly,
    "file"=FaFile(assembly),
    "bsgenome"=assembly
  )

  return(assembly)

}
