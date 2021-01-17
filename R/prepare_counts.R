#' Format Counts
#'
#' @description
#' Format TSS or TSR counts for further analysis. 
#'
#' @inheritParams common_params
#' @param data_type Whether to format TSS or TSR counts 
#'
#' @details
#' When TSSs or TSRs are first loaded into the TSRexploreR object
#'   they are stored as GRanges objects.
#' This function converts these into data.table format,
#'   and adds a few important columns for downstream analysis.
#'
#' @return TSRexploreR object with properly formatted features
#'    in data.table format.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#'
#' @rdname format_counts-function
#' @export

format_counts <- function(
  experiment,
  data_type=c("tss", "tsr"),
  samples="all"
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr"))
  assert_that(is.character(samples))

  ## Get selected samples and generate raw count matrices.
  if (data_type == "tss") {

    ## Get selected samples.
    select_samples <- tss_experiment(experiment)
    if (!all(samples == "all")) select_samples <- select_samples[samples]

  } else if (data_type == "tsr") {

    ## Get selected samples.
    select_samples <- tsr_experiment(experiment)
    if (!all(samples == "all")) select_samples <- select_samples[samples]

  }

  ## Turn counts into data.table
  if (data_type %in% c("tss", "tsr")) {
    raw_counts <- map(select_samples, function(x) {
      x <- as.data.table(x)
      x[,
        FHASH := str_c(seqnames, start, end, strand, sep=":"),
        by=seq_len(nrow(x))
      ]
      return(x)
    })
  }

  ## Place counts in proper TSRexploreR object slot.
  if (data_type == "tss") {
    experiment@counts$TSSs$raw <- c(experiment@counts$TSSs$raw, raw_counts)
  } else if (data_type == "tsr") {
    experiment@counts$TSRs$raw <- c(experiment@counts$TSRs$raw, raw_counts)
  }

  return(experiment)
}

#' Feature Counts
#'
#' Count the number of reads associated with genes or transcripts
#'   based on the aggregate score of TSSs or TSRs annotated to them.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'
#'
#' @details
#' The 'count_features' function counts the total score of TSSs or TSRs associated
#'   with a gene or transcript.
#' This allows for an RNA-seq like analysis of gene expression using TSS mapping data.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data=annotation,
#'   data_type="tss", feature_type="transcript"
#' )
#' tsre_exp <- count_features(tsre_exp, data_type="tss")
#'
#' @return RNA-seq like count matrix added to tsr explorer object.
#'
#' @seealso
#' \code{\link{annotate_features}} to first annotate the TSSs or TSRs.
#'
#' @rdname count_features-function
#' @export

count_features <- function(
  experiment,
  data_type=c("tss", "tsr")
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))  
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr"))
  
  ## Get information on whether annotation was by gene or transcript.
  anno_type <- ifelse(
    experiment@settings$annotation[, feature_type] == "transcript",
    "transcriptId", "geneId"
  )

  ## Extract appropriate counts.
  sample_data <- extract_counts(experiment, data_type, "all")

  ## Get feature counts.
  sample_data <- sample_data %>%
    map(function(x) {
      feature_scores <- x[
        !simple_annotations %in% c("Downstream", "Intergenic", "Antisense"),
        .(feature_score=sum(score)),
        by=eval(anno_type)
      ]

      setkeyv(x, anno_type)
      setkeyv(feature_scores, anno_type)
      x <- merge(x, feature_scores, all.x=TRUE)
      setcolorder(x, c(discard(colnames(x), ~ . == anno_type), anno_type))
      
      return(x)
    })

  ## Get feature counts per sample.
  counts <- sample_data %>%
    map(function(x) {
      setnames(x, old=anno_type, new="feature")
      x <- unique(x[, .(feature, feature_score)])
      x[, feature_score := ifelse(is.na(feature_score), 0, feature_score)]
      return(x)
    })

  ## Store counts in appropriate slots.
  walk(counts, ~ setnames(., old=c("feature", "feature_score"), new=c(anno_type, "score")))
  walk(sample_data, ~ setnames(., old="feature", new=anno_type))

  if (data_type == "tss") {
    experiment@counts$TSSs$raw <- sample_data
    experiment@counts$TSS_features$raw <- counts
  } else {
    experiment@counts$TSRs$raw <- sample_data
    experiment@counts$TSR_features$raw <- exp_counts
  }

  return(experiment)
}

#' Count Matrix
#'
#' Generate count matrices for correlation analysis.
#'
#' @param experiment TSRexploreR object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param samples Character vector of samples to turn into a count matrix
#'
#' @details
#' In order to compare samples, a matrix must first be generated with
#'   genomic positions as rows and samples as columns.
#'
#' For each TSS, the genomic position is simply its single base position
#' However, for TSRs, overlapping TSRs between samples are first merged
#'   to generate consensus TSRs for all considered samples.
#'
#' @return TSRexploreR object with added count matrices
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' tsre_exp <- count_matrix(tsre_exp, data_type="tss")
#'
#' @seealso \code{\link{tmm_normalize}} to normalize the count matrices.
#'   \code{\link{plot_correlation}} for various correlation plots.
#'
#' @rdname count_matrix-function

count_matrix <- function(
  experiment,
  data_type=c("tss", "tsr", "tss_features", "tsr_features"),
  samples="all"
) {
  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features")
  )
  assert_that(is.character(samples)) 

  ## Extract counts.
  select_samples <- extract_counts(experiment, data_type, samples)

  if (data_type %in% c("tss", "tsr")) {
    select_samples <- map(select_samples, function(x) {
      x <- x[, .(seqnames, start, end, strand, score, FHASH)]
      return(x)
    })
  }

  ## Generate count matrix.
  if (data_type == "tss") {

    ## Merge overlapping TSSs.
    select_samples <- rbindlist(select_samples, idcol="sample")
    select_samples <- dcast(
      select_samples,
      seqnames + start + end + strand + FHASH ~ sample,
      value.var="score", fill=0
    )

  } else if (data_type == "tsr") {

    ## Merge overlapping TSRs.
    tsr_consensus <- select_samples %>%
      map(as_granges, keep_mcols=FALSE) %>%
      bind_ranges %>%
      reduce_ranges_directed %>%
      as.data.table(key=c("seqnames", "strand", "start", "end"))

    ## Create raw count matrix.
    select_samples <- select_samples %>%
      imap(function(x, y) {
        consensus_counts <- x %>%
          foverlaps(tsr_consensus) %>%
          {.[, .(score=sum(score)), by=.(seqnames, start, end, strand)]}
        setnames(consensus_counts, old="score", new=y)
      }) %>%
      purrr::reduce(
        merge, by=c("seqnames", "start", "end", "strand"),
        all=TRUE
      )
    select_samples[is.na(select_samples)] <- 0

    select_samples[,
      FHASH := str_c(seqnames, start, end, strand, sep=":"),
      by=seq_len(nrow(select_samples))
    ]


  } else if (data_type %in% c("tss_features", "tsr_features")) {
    
    ## Get annotation type.
    anno_type <- experiment@settings$annotation[["feature_type"]]
    anno_type <- ifelse(anno_type == "gene", "geneId", "transcriptId")

    ## Change feature counts to matrix.
    select_samples <- rbindlist(select_samples, idcol="sample")
    select_samples <- dcast(
      select_samples,
      as.formula(str_c(anno_type, " ~ sample")),
      fill=0
    )

  }

  ## Create SummarizedExperiments.
  if (data_type %in% c("tss", "tsr")) {
    
    ## Prepare data for RangedSummarizedExperiment
    row_ranges <- as_granges(select_samples, keep_mcols=FALSE)

    select_samples[, c("seqnames", "start", "end", "strand") := NULL]
    select_samples <- as.matrix(select_samples, rownames="FHASH")

    col_data <- DataFrame(
      samples=colnames(select_samples),
      row.names=colnames(select_samples)
    )
  
    ## Make RangedSummarizedExperiment.
    rse <- SummarizedExperiment(
      assays=list(counts=select_samples),
      rowRanges=row_ranges,
      colData=col_data
    )
  } else if (data_type %in% c("tss_features", "tsr_features")) {
    
     ## Prepare data for SummarizedExperiment.
    row_data <- select_samples[, 1]
    
    select_samples[, eval(anno_type) := NULL]
    select_samples <- as.matrix(select_samples)
    rownames(select_samples) <- row_data[[1]]

    col_data <- DataFrame(
      samples=colnames(select_samples),
      row.names=colnames(select_samples)
    )

    ## Make SummarizedExperiment.
    rse <- SummarizedExperiment(
      assays=list(counts=select_samples),
      rowData=row_data,
      colData=col_data
    )
  }

  ## Add RangedSummarizedExperiment back to trexplorer object.
  experiment <- set_count_slot(
    experiment, rse,
    "counts", data_type, "matrix"
  )

  return(experiment)  
}
