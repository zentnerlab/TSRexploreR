#' Genes/Transcripts Detected
#'
#' @description
#' This plotting function returns a stacked barplot showing the number of
#' features detected with and without a promoter proximal TSS or TSR. The information 
#'
#' @inheritParams common_params
#' @param data_type Whether TSSs or TSRs should be analyzed.
#' @param ... Arguments passed to geom_col.
#'
#' @details
#' This function will returnthe number of genes or transcripts with an associated 
#' unique TSS or TSR. Information on whether the feature has a promoter-proximal 
#' TSS or TSR is included in the output for plotting purposes.
#'
#' A set of functions to control data structure for plotting are included. 'use_normalized' 
#' will use  normalized scores, which only matters if 'consider_score' is TRUE.
#' 'threshold' defines the minimum number of raw counts a TSS or TSR must have to be 
#' considered. dominant' specifies whether only the dominant TSS or TSR (determined
#' using the 'mark_dominant' function) is considered. For TSSs, this can be either 
#' dominant TSS per TSR or gene/transcript, and for TSRs it is the dominant TSR 
#' per gene/transcript. 'data_conditionals' can be used to filter, quantile, order, 
#' and/or group data for plotting.
#'
#' @return ggplot2 object of detected features.
#'
#' @seealso
#' \code{\link{annotate_features}} to annotate the TSSs or TSRs.
#'
#' @examples
#' library("magrittr")
#' data(TSSs)
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer(genome_annotation=annotation) %>%
#'   format_counts(data_type="tss") %>%
#'   annotate_features(data_type="tss")
#'
#' # Detected features with TSSs plot.
#' \donttest{plot_detected_features(tsre, data_type="tss")}
#'
#' # Detected features with TSRs plot.
#' tsre <- tsre %>%
#'   tss_clustering(threshold=3) %>%
#'   annotate_features(data_type="tss")
#' \donttest{plot_detected_features(tsre, data_type="tsr")}
#'
#' @export

plot_detected_features <- function(
  experiment,
  samples="all",
  data_type=c("tss", "tsr"),
  threshold=NULL,
  dominant=FALSE,
  use_normalized=FALSE,
  data_conditions=NULL,
  return_table=FALSE,
  ...
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.flag(dominant))
  assert_that(is.flag(use_normalized))
  assert_that(is.null(data_conditions) || is.list(data_conditions))
  assert_that(is.flag(return_table))

  ## Get sample data.
  sample_data <- experiment %>%
    extract_counts(data_type, samples, use_normalized) %>%
    preliminary_filter(dominant, threshold)
  
  ## Apply data conditioning if requested.
  sample_data <- condition_data(sample_data, data_conditions)

  ## Rename feature column.
  walk(sample_data, function(x) {
    setnames(
      x, old=ifelse(
        experiment@settings$annotation[, feature_type] == "transcript",
        "transcriptId", "geneId"
      ),
      new="feature"
    )
  })

  ## Get feature counts.
  grouping_status <- case_when(
    !is.null(data_conditions$quantiling) ~ "row_quantile",
    !is.null(data_conditions$grouping) ~ "row_groups",
    TRUE ~ "none"
  )

  sample_data <- rbindlist(sample_data, idcol="sample")
  sample_data <- .count_features(sample_data, grouping_status)

  ## Prepare data for plotting.
  sample_data[, total := NULL]
  plot_data <- melt(
    sample_data,
    measure.vars=c("with_promoter", "without_promoter"),
    variable.name="count_type",
    value.name="feature_count"
  )
  plot_data[, count_type := factor(
    count_type, levels=c("without_promoter", "with_promoter")
  )]

  ## Order samples if required.
  if (!all(samples == "all")) {
    plot_data[, sample := factor(sample, levels=samples)]
  }

  ## Return a table if required.
  if (return_table) return(as.data.frame(plot_data))

  ## Plot data.
  p <- ggplot(plot_data, aes(x=.data$sample, y=.data$feature_count, fill=.data$count_type)) +
    geom_col(position="stack", ...) +
    theme_bw() +
    ylim(c(0, NA)) +
    ylab("Feature Count") +
    xlab("Sample") +
    theme(
      axis.text.x=element_text(angle=45, hjust=1)
    )

  if (grouping_status != "none") {
    p <- p + facet_grid(fct_rev(factor(grouping)) ~ .)
  }

  return(p)

}

#' Calculate Feature Counts
#'
#' @param sample_data Sample data.
#' @param grouping_status Whether there is a grouping variable.

.count_features <- function(sample_data, grouping_status) {

  if (grouping_status != "none") {
    setnames(sample_data, old=grouping_status, new="grouping")
    sample_data <- sample_data[,
      .(grouping, promoter=any(simple_annotations == "Promoter")),
      by=.(sample, feature)
    ][,
      .(with_promoter=sum(promoter), without_promoter=.N - sum(promoter), total=.N),
      by=.(sample, grouping)
    ]
  } else {
    sample_data <- sample_data[,
      .(promoter=any(simple_annotations == "Promoter")),
      by=.(sample, feature)
    ][,
      .(with_promoter=sum(promoter), without_promoter=.N - sum(promoter), total=.N),
      by=sample
    ]
  }

  return(sample_data)
}
