#' Mark Dominant
#'
#' @description
#' Mark TSSs as dominant per TSR or gene/transcript, or TSRs as dominant per gene/transcript. 
#'
#' @inheritParams common_params
#' @param data_type Whether to mark dominant TSSs or TSRs.
#' @param mark_per By default, the function marks dominant TSR per gene, and dominant TSS per TSR.
#'   TSSs can also be set per as dominant TSS per 'gene'.
#'
#' @details
#' This function marks which TSSs are dominant per TSR or gene, or which TSR is
#' dominant per gene. Analysis of dominant features may help to cut through the 
#' noise to get information such as the length of the primary 5' UTR and sequence
#' features associated with the the strongest TSS.
#'
#' Setting a 'threshold' will only mark a TSS or TSR as dominant if its score 
#' is greater than or equal to the threshold.
#'
#' 'mark_per' controls the behavior of the function. For TSSs 'default' will mark,
#' dominant TSS per TSR, and for TSRs the dominant TSR per gene is marked. For TSSs, 
#' 'gene' can also be specified, which will mark the dominant TSS per gene.  
#'
#' @return TSRexploreR object with dominant status added to TSSs or TSRs.
#'
#' @seealso
#' \code{\link{associate_with_tsr}} to associate TSSs with TSRs prior to marking
#'   dominant TSS per TSR.
#'
#' @examples
#' data(TSSs)
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer(genome_annotation=annotation) %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3) %>%
#'   associate_with_tsr
#'
#' # Dominant TSS per gene.
#' tsre <- annotate_features(tsre, data_type="tss")
#' \donttest{mark_dominant(tsre, data_type="tss", mark_per="gene")}
#'
#' # Dominant TSS per TSR.
#' \donttest{mark_dominant(tsre, data_type="tss")}
#'
#' # Dominant TSR per gene.
#' tsre <- annotate_features(tsre, data_type="tsr")
#' \donttest{mark_dominant(tsre, data_type="tsr")}
#'
#' @export

mark_dominant <- function(
  experiment,
  data_type=c("tss", "tsr"),
  threshold=NULL,
  use_normalized=FALSE,
  mark_per="default",
  exclude_antisense=TRUE
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr"))
  assert_that(
    is.null(threshold) ||
    (is.numeric(threshold) && threshold >= 0)
  )
  mark_per <- match.arg(str_to_lower(mark_per), c("default", "gene"))

  ## Select samples.
  select_samples <- extract_counts(experiment, data_type, "all", use_normalized)

  ## Set threshold to 0 if not supplied.
  if (is.null(threshold)) threshold <- 0

  ## Mark whether antisense should be excluded
  anno_remove <- ifelse(
    exclude_antisense,
    c("Downstream", "Intergenic", "Antisense"),
    c("Downstream", "Intergenic")
  )

  ## Mark dominant TSS per gene if requested.
  if (data_type == "tss" & mark_per == "gene") {
    dominant <- map(select_samples, function(x) {
      x[
        !is.na(TSR_FHASH),
        dominant := (
          score == max(score) &
          !simple_annotations %in% anno_remove &
          score >= threshold
        ),
        by=eval(ifelse(
          experiment@settings$annotation[, feature_type] == "transcript",
          "transcriptId", "geneId"
        ))
      ]

      return(x)
    })
    
  ## Mark the dominant TSR per gene if requested.
  } else if (data_type == "tsr") {
    dominant <- map(select_samples, function(x) {
      x[,
        dominant := (
          score == max(score) &
          !simple_annotations %in% anno_remove &
          score >= threshold
        ),
        by=eval(ifelse(
          experiment@settings$annotation[, feature_type] == "transcript",
          "transcriptId", "geneId"
        ))
      ]

      return(x)
    })
    
  ## Mark the dominant TSS per TSR if requested.
  } else if (data_type == "tss" & mark_per == "default") {
    dominant <- map(select_samples, function(x) {
      x[
        !is.na(TSR_FHASH),
        dominant := (
          !is.na(score) &
          score == max(score) &
          score >= threshold
        ),
        by=TSR_FHASH
      ]

      return(x)
    })
  }

  ## Return dominant TSS/TSR.
  experiment <- set_count_slot(
    experiment, dominant,
    "counts", data_type, "raw"
  )

  return(experiment)
}
