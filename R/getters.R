
#' Get TSS or TSR GRanges.
#'
#' @description
#' Extract the TSS or TSR GRanges from the 'experiment' slot.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#'
#' @return List of TSS or TSR GRanges.
#'
#' @seealso
#' \code{\link{get_counts}} to get TSS or TSR counts data.frame.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#'
#' tsre <- tsr_explorer(TSSs[1])
#' get_granges(tsre)
#'
#' @export

get_granges <- function(
  experiment,
  data_type=c("tss", "tsr"),
  samples="all"
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr")
  )
  assert_that(is.character(samples))

  ## Get sample names.
  if (samples == "all") {
    samples <- switch(
      data_type,
      "tss"=names(experiment@experiment$TSSs),
      "tsr"=names(experiment@experiment$TSRs)
    )
  }

  ## Retrieve GRanges.
  grs <- switch(
    data_type,
    "tss"=experiment@experiment$TSSs,
    "tsr"=experiment@experiment$TSRs
  )

  return(grs)
}

#' Get TSS or TSR count tables.
#'
#' @description
#' Extract the TSS or TSR count tables from the 'counts' slot.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#'
#' @return List of TSS or TSR count data.frame.
#'
#' @seealso
#' \code{\link{get_granges}} to get TSS or TSR GRanges.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#'
#' tsre <- tsr_explorer(TSSs[1]) %>%
#'   format_counts(data_type="tss")
#' get_counts(tsre)
#'
#' @export

get_counts <- function(
  experiment,
  data_type=c("tss", "tsr"),
  samples="all"
) {
  
  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr")
  )
  assert_that(is.character(samples))

  ## Get sample names.
  if (samples == "all") {
    samples <- switch(
      data_type,
      "tss"=names(experiment@counts$TSSs$raw),
      "tsr"=names(experiment@counts$TSRs$raw)
    )
  }
  
  ## Retrieve GRanges.
  cts <- switch(
    data_type,
    "tss"=experiment@counts$TSSs$raw,
    "tsr"=experiment@counts$TSRs$raw
  )
  cts <- map(cts, as.data.frame)

  return(cts)
}

#' Get Genome Annotation
#'
#' Extract the genome annotation.
#'
#' @inheritParams common_params
#'
#' @return TxDb object of genome annotation.
#'
#' @seealso
#' \code{\link{get_assembly}} to get genome assembly.
#'
#' @examples
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#' tsre <- tsr_explorer(genome_annotation=annotation)
#' get_annotation(tsre)
#'
#' @export

get_annotation <- function(
  experiment
) {

  ## Input check.
  assert_that(is(experiment, "tsr_explorer"))

  ## Get genome annotation.
  an <- experiment@meta_data$genome_annotation

  return(an)
}

#' Get Genome Assembly
#'
#' Extract the genome assembly.
#'
#' @inheritParams common_params
#'
#' @return Either BSgenome or FaFile depending on original source of assembly.
#'
#' @seealso
#' \code{\link{get_annotation}} to get genome annotation.
#'
#' @examples
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#' tsre <- tsr_explorer(genome_assembly=assembly)
#' get_assembly(tsre)
#'
#' @export

get_assembly <- function(
  experiment
) {

  ## Input check.
  assert_that(is(experiment, "tsr_explorer"))

  ## Get genome annotation.
  as <- experiment@meta_data$genome_assembly

  return(as)
}

#' Get Sample Sheet
#'
#' Extract the sample sheet.
#'
#' @inheritParams common_params
#'
#' @return data.frame
#'
#' @examples
#' sample_sheet <- data.frame(
#'   sample_name=sprintf("S288C_WT_%s", seq_len(3)),
#'   file_1=NA, file_2=NA,
#'   condition=rep("Untreated", 3)
#' )
#' tsre <- tsr_explorer(sample_sheet=sample_sheet)
#' get_sample_sheet(tsre)
#'
#' @export

get_sample_sheet <- function(
  experiment
) {

  ## Input check.
  assert_that(is(experiment, "tsr_explorer"))

  ## Get the sample sheet.
  ss <- experiment@meta_data$sample_sheet

  return(as.data.frame(sample_sheet))
}

#' Get Shifting Results
#'
#' Extract the TSS cluster shifting results.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#'
#' @return List of TSS cluster shifting data.frames.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "TSRexploreR")
#' TSSs <- readRDS(TSSs)
#'
#' tsre <- TSSs[c(1, 4)] %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3) %>%
#'   tss_shift(
#'     tsre,
#'     sample_1=c(TSS="S288C_WT_1", TSR="S288C_WT_1"),
#'     sample_2=c(TSS="S288C_D_1", TSR="S288C_D_1"),
#'     comparison_name="Untreated_vs_Diamide",
#'     max_distance = 100, min_threshold = 10, n_resamples = 1000L
#'   )
#' get_shifting_results(tsre)
#'
#' @export

get_shifting_results <- function(
  experiment,
  samples="all"
) {

  ## Input check.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))

  ## Get sample names.
  if (samples == "all") {
    samples <- names(samples@shifting$results)
  }

  ## Get shifting results.
  shft <- experiment@shifting$results[samples]
  shft <- map(shft, as.data.frame)

  return(shft)
}

#' Get Differential TSS/TSR models
#'
#' Get the DESeq2 or edgeR model.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#'
#' @return DESeq2 or edgeR differential expression model.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=NA, file_2=NA,
#'   condition=c(
#'     rep("Diamide", 3),
#'     rep("Untreated", 3)
#'   )
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet) %>%
#'   format_counts(data_type="tss") %>%
#'   fit_de_model(tsre, ~condition, data_type="tss")
#' get_diff_model(tsre, data_type="tss")
#'
#' @export

get_diff_model <- function(
  experiment,
  data_type=c("tss", "tsr")
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr")
  )

  ## Get diff expression model.
  dem <- switch(
    data_type,
    "tss"=experiment@diff_features$TSSs$model,
    "tsr"=experiment@diff_features$TSRs$model
  )

  return(dem)
}

#' Get Differential TSS/TSR results
#'
#' Get the DESeq2 or edgeR differential features.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#'
#' @return List of differential features data.frames.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=NA, file_2=NA,
#'   condition=c(
#'     rep("Diamide", 3),
#'     rep("Untreated", 3)
#'   )
#' )
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet) %>%
#'   format_counts(data_type="tss") %>%
#'   fit_de_model(~condition, data_type="tss") %>%
#'   differential_expression(
#'     tsre, data_type="tss",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="contrast",
#'     comparison=c("condition", "Diamide", "Untreated")
#'   )
#' get_diff_results(tsre, data_type="tss")
#'
#' @export

get_diff_results <- function(
  experiment,
  data_type=c("tss", "tsr"),
  samples="all"
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr")
  )
  assert_that(is.character(samples))

  ## Get the sample names.
  if (samples == "all") {
    samples <- switch(
      data_type,
      "tss"=names(experiment@diff_features$TSSs$results),
      "tsr"=names(experiment@diff_features$TSRs$results)
    )
  }

  ## Get the differential expression results.
  der <- switch(
    data_type,
    "tss"=experiment@diff_features$TSSs$results[samples],
    "tsr"=experiment@diff_features$TSRs$results[samples]
  )
  der <- map(der, as.data.frame)

  return(der)
}