
#' Get TSS or TSR GRanges.
#'
#' @description
#' Extract the TSS or TSR GRanges from the 'experiment' slot.
#'
#' @inheritParams common_params
#' @param data_type Get either TSS ('tss') or TSR ('tsr') GRanges.
#'
#' @return List of TSS or TSR GRanges.
#'
#' @seealso
#' \code{\link{get_counts}} to get TSS or TSR counts data.frame.
#'
#' @examples
#' data(TSSs_reduced)
#'
#' exp <- tsr_explorer(TSSs_reduced)
#' 
#' gr <- get_granges(exp)
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
#' @param data_type Get either TSS ('tss') or TSR ('tsr') count tables.
#'
#' @return List of TSS or TSR count data.frame.
#'
#' @seealso
#' \code{\link{get_granges}} to get TSS or TSR GRanges.
#'
#' @examples
#' data(TSSs_reduced)
#'
#' exp <- tsr_explorer(TSSs_reduced) %>%
#'   format_counts(data_type="tss")
#' 
#' cts <- get_counts(exp)
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
#' exp <- tsr_explorer(genome_annotation=annotation)
#' 
#' get_annotation(exp)
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
#' exp <- tsr_explorer(genome_assembly=assembly)
#' 
#' a <- get_assembly(exp)
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
#' 
#' exp <- tsr_explorer(sample_sheet=sample_sheet)
#' 
#' ss <- get_sample_sheet(exp)
#'
#' @export

get_sample_sheet <- function(
  experiment
) {

  ## Input check.
  assert_that(is(experiment, "tsr_explorer"))

  ## Get the sample sheet.
  ss <- experiment@meta_data$sample_sheet

  return(as.data.frame(ss))
}

#' Get Shifting Results
#'
#' Extract the TSR shifting results.
#'
#' @inheritParams common_params
#'
#' @return List of TSR shifting data.frames.
#'
#' @examples
#' data(TSSs)
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "TSRexploreR")
#' samples <- data.frame(
#'   sample_name=c(sprintf("S288C_D_%s", seq_len(3)), sprintf("S288C_WT_%s", seq_len(3))),
#'   file_1=rep(NA, 6), file_2=rep(NA, 6),
#'   condition=c(rep("Diamide", 3), rep("Untreated", 3))
#' )
#'
#' exp <- TSSs %>%
#'   tsr_explorer(sample_sheet=samples, genome_assembly=assembly) %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3) %>%
#'   merge_samples(data_type = "tss", merge_group="condition") %>%
#'   merge_samples(data_type = "tsr", merge_group="condition") %>%
#'   tss_shift(
#'     sample_1=c(TSS="S288C_WT_1", TSR="S288C_WT_1"),
#'     sample_2=c(TSS="S288C_D_1", TSR="S288C_D_1"),
#'     comparison_name="Untreated_vs_Diamide",
#'     max_distance = 100, min_threshold = 10, n_resamples = 1000L
#'   )
#'   
#' sr <- get_shifting_results(exp)
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
  if (all(samples == "all")) {
    samples <- names(experiment@shifting$results)
  }

  ## Get shifting results.
  shft <- experiment@shifting$results[samples]

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
#' data(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=rep(NA, 6), file_2=rep(NA, 6),
#'   condition=c(
#'     rep("Diamide", 3),
#'     rep("Untreated", 3)
#'   )
#' )
#'
#' exp <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet) %>%
#'   format_counts(data_type="tss") %>%
#'   fit_de_model(~condition, data_type="tss", method="edgeR")
#' 
#' dm <- get_diff_model(exp, data_type="tss")
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
#' data(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=c(
#'     sprintf("S288C_D_%s", seq_len(3)),
#'     sprintf("S288C_WT_%s", seq_len(3))
#'   ),
#'   file_1=rep(NA, 6), file_2=rep(NA, 6),
#'   condition=c(
#'     rep("Diamide", 3),
#'     rep("Untreated", 3)
#'   )
#' )
#'
#' exp <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet) %>%
#'   format_counts(data_type="tss") %>%
#'   fit_de_model(~condition, data_type="tss", method="edgeR") %>%
#'   differential_expression(
#'     data_type="tss",
#'     comparison_name="Diamide_vs_Untreated",
#'     comparison_type="coef",
#'     comparison=2
#'   )
#' 
#' dr <- get_diff_results(exp, data_type="tss")
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
  if (all(samples == "all")) {
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
