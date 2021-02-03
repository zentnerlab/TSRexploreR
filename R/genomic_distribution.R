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
#' This plotting function will create a stacked barplot of genomic feature types for each sample.
#' Genomic features include exons, introns, intergenic, downstream, antisense, and promoter regions.
#' The promoter region is user-defined during annotation.
#'
#' A set of functions to control data structure for plotting are included. 'use_normalized' 
#' will use  normalized scores, which only matters if 'consider_score' is TRUE.
#' 'threshold' defines the minimum number of raw counts a TSS or TSR must have to be 
#' considered. dominant' specifies whether only the dominant TSS or TSR (determined
#' using the 'mark_dominant' function) is considered. For TSSs, this can be either 
#' dominant TSS per TSR or gene/transcript, and for TSRs it is the dominant TSR 
#' per gene/transcript. 'data_conditions' can be used to filter, quantile, order, 
#' and/or group data for plotting.
#'
#' If 'return_table' is TRUE, a data.frame containing the underlying data
#'   for the plot is returned.
#'
#' @return ggplot2 plot with TSS or TSR genomic distribution.
#'   If 'return_table' is TRUE returns a data.frame of underlying stats.
#'
#' @seealso \code{\link{annotate_features}} to annotate TSSs or TSRs.
#'
#' @examples
#' data(TSSs)
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer(genome_annotation=annotation) %>%
#'   format_counts(data_type="tss") %>%
#'   annotate_features(data_type="tss")
#'
#' # TSS genomic distribution plot.
#' \donttest{plot_genomic_distribution(tsre, data_type="tss")}
#'
#' # TSR genomic distribution plot.
#' tsre <- tsre %>%
#'   tss_clustering(threshold=3) %>%
#'   annotate_features(data_type="tsr")
#' \donttest{plot_genomic_distribution(tsre, data_type="tsr")}
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
    genomic_dist[, sample := factor(sample, levels=rev(samples))]
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
