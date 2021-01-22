#' Genomic Distribution
#'
#' @description
#' Get genomic distribution of TSSs and TSRs.
#'
#' @inheritParams common_params
#' @param data_type Whether to get distribution of TSSs or TSRs.
#' @param ... Arguments passed to geom_col.
#'
#' @details
#' This function summarizes the distribution of TSSs or TSRs relative to annotated 
#' genomic features (exons, introns, intergenic, downstream, and promoter regions).
#' The promoter region is user-defined during annotation.
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
#' @return DataFrame with TSS or TSR genomic distribution stats.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#' exp <- format_counts(exp, data_type="tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#' exp <- annotate_features(
#'   exp, annotation_data=annotation,
#'   data_type="tss", feature_type="transcript"
#' )
#' genomic_dist <- genomic_distribution(exp, data_type="tss")
#'
#' @seealso \code{\link{annotate_features}} to annotate TSSs or TSRs.
#'   \code{\link{plot_genomic_distribution}} to plot the genomic distribution.
#'
#' @export 

plot_genomic_distribution <- function(
  experiment,
  data_type=c("tss", "tsr", "shift"),
  samples="all",
  threshold=NULL,
  use_normalized=FALSE,
  dominant=FALSE,
  data_conditions=NULL,
  return_table=FALSE,
  ...
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "shift"))
  assert_that(is.character(samples))
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.flag(dominant))
  assert_that(is.null(data_conditions) || is.list(data_conditions))
  assert_that(is.flag(return_table))

  ## Get samples.
  selected_samples <- experiment %>%
    extract_counts(data_type, samples, use_normalized) %>%
    preliminary_filter(dominant, threshold)

  walk(selected_samples, function(x) {
    x[, simple_annotations := factor(
      simple_annotations,
      levels=c("Promoter", "Exon", "Intron", "Downstream", "Intergenic", "Antisense")
    )]
  })

  ## Apply advanced grouping.
  selected_samples <- condition_data(selected_samples, data_conditions)

  ## Calculate distribution.
  selected_samples <- rbindlist(selected_samples, idcol="sample")

  grouping_status <- case_when(
    !is.null(data_conditions$quantiling) ~ "row_quantile",
    !is.null(data_conditions$grouping) ~ "row_groups",
    TRUE ~ "none"
  )

  genomic_dist <- .calculate_distribution(selected_samples, grouping_status)

  ## Order samples if required.
  if (!all(samples == "all")) {
    genomic_dist[, sample := factor(sample, levels=samples)]
  }

  ## Return table is requested.
  if (return_table) return(as.data.frame(genomic_dist))

  ## Plot the genomic distribution.
  p <- ggplot(genomic_dist, aes(x=.data$sample, y=.data$count, fill=fct_rev(.data$simple_annotations))) +
    geom_col(position="fill", ...) +
    coord_flip() +
    ylab("Fraction") +
    theme_bw() +
    theme(
      axis.title.y=element_blank(),
      panel.grid=element_blank()
    )

  if (grouping_status != "none") {
    p <- p + facet_grid(grouping ~ .)
  }

  return(p)

}

#' Calculate Genomic Distribution
#'
#' @param selected_samples Samples to analyze.
#' @param grouping_status Whether data is grouped.

.calculate_distribution <- function(selected_samples, grouping_status) {

  if (grouping_status != "none") {
    setnames(selected_samples, old=grouping_status, new="grouping")
    genomic_distribution <- selected_samples[,
      .(count=.N),
      by=.(sample, simple_annotations, grouping)
    ][,
      .(simple_annotations, count, fraction=count / sum(count)),
      by=.(sample, grouping)
    ]
  } else {
    genomic_distribution <- selected_samples[,
      .(count=.N),
      by=.(sample, simple_annotations)
    ][,
      .(simple_annotations, count, fraction=count / sum(count)),
      by=sample
    ]
  }

  return(genomic_distribution)

}

#' Plot Genomic Distribution
#'
#' Plot genomic distribution of TSSs or TSRs.
#'
#' @param genomic_distribution Dataframe of TSS or TSR genomic distributions from tsr_genomic_distribution.
#'
#' @details
#' This plotting function will create a stacked barplot of genomic feature types for each sample.
#' The underlying data used for plotting comes from the 'genomic_distribution' function.
#' Genomic features include exons, introns, intergenic, downstream and promoter regions.
#' The promoter region is user-defined during annotation.
#'
#' @return ggplot2 object with TSS or TSR genomic distribution plot.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#' exp <- format_counts(exp, data_type="tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#' exp <- annotate_features(
#'   exp, annotation_data=annotation,
#'   data_type="tss", feature_type="transcript"
#' )
#' genomic_dist <- genomic_distribution(exp, data_type="tss")
#' plot_genomic_distribution(genomic_dist)
#'
#' @seealso \code{\link{annotate_features}} to annotate TSSs or TSRs.
#'   \code{\link{genomic_distribution}} to prepare annotated TSSs or TSRs for plotting.

plot_genomic_dist <- function(x) NULL
