#' Genomic Distribution
#'
#' @description
#' Get genomic distribution of TSSs and TSRs
#'
#' @param experiment TSRexploreR object with annotated TSRs
#' @param samples Either "all" or a vector of sample names
#' @param data_type Whether to get distribution of TSSs or TSRs
#' @param threshold Raw count threshold for a TSS or TSR to be considered
#' @param dominant Whether only the dominant TSS per gene or TSR should be considered
#' @param data_conditions Apply conditions to data (supports filtering and quantiles/grouping)
#' @param use_normalized Whether normalized or raw counts should be used
#'
#' @details
#' This function summarizes the distribution of TSSs or TSRs relative
#'   to annotated genomic features (exons, introns, intergenic, 
#'   downstream and promoter regions.)
#' The promoter region is user defined when annotating the TSSs or TSRs.
#'
#' A set of functions to control data structure for plotting are included.
#' 'threshold' defines the minimum number of raw counts a TSS or TSR
#'  must have to be considered.
#' 'dominant' specifies whether only the dominant TSS or TSR is considered 
#'   (annotated by the 'mark_dominant' function).
#' For TSSs this can be either dominant per TSR or gene, and for TSRs
#'   it is just the dominant TSR per gene.
#' 'data_conditions' allows for the advanced filtering, ordering, and grouping
#'   of data.
#'
#' @return DataFrame with TSS or TSR genomic distribution stats
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
#' genomic_dist <- genomic_distribution(tsre_exp, data_type="tss")
#'
#' @seealso \code{\link{annotate_features}} to annotate TSSs or TSRs.
#'   \code{\link{plot_genomic_distribution}} to plot the genomic distribution.
#'
#' @rdname genomic_distribution-function
#' @export 

genomic_distribution <- function(
  experiment,
  data_type=c("tss", "tsr"),
  samples="all",
  threshold=NULL,
  use_normalized=FALSE,
  dominant=FALSE,
  data_conditions=NA
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr"))
  assert_that(is.character(samples))
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.flag(dominant))
  if (all(!is.na(data_conditions)) && !is(data_conditions, "list")) stop("data_conditions must in list form")

  ## Extract samples.
  selected_samples <- experiment %>%
    extract_counts(data_type, samples, use_normalized) %>%
    preliminary_filter(dominant, threshold)

  walk(selected_samples, function(x) {
    x[, simple_annotations := factor(
      simple_annotations,
      levels=c("Promoter", "Exon", "Intron", "Downstream", "Intergenic")
    )]
  })

  ## Apply advanced grouping.
  if (all(!is.na(data_conditions))) {
    selected_samples <- do.call(group_data, c(list(signal_data=selected_samples), data_conditions))
  }

  ## Prepare data to be plotted later.
  selected_samples <- rbindlist(selected_samples, idcol="samples")
  groupings <- any(names(data_conditions) %in% c("quantile_by", "grouping"))
  genmic_distribution <- .calculate_distribution(selected_samples, groupings)

  ## Order samples if required.
  if (!all(samples == "all")) {
    genomic_distribution[, samples := factor(samples, levels=samples)]
  }

  ## Prepare dataframe to return.
  dist_exp <- DataFrame(genomic_distribution)

  ## Add quantile information to summarized experiment.
  metadata(dist_exp)$groupings <- groupings
  metadata(dist_exp)$dominant <- dominant

  return(dist_exp)
}

#' Calculate Genomic Distribution
#'
#' @param selected_samples Samples to analyze
#' @param groupings Whether data is grouped

.calculate_distribution <- function(selected_samples, groupings) {

  if (groupings) {
    genomic_distribution <- selected_samples[,
      .(count=.N),
      by=.(samples, simple_annotations, grouping)
    ][,
      .(simple_annotations, count, fraction=count / sum(count)),
      by=.(samples, grouping)
    ]
  } else {
    genomic_distribution <- selected_samples[,
      .(count=.N),
      by=.(samples, simple_annotations)
    ][,
      .(simple_annotations, count, fraction=count / sum(count)),
      by=samples
    ]
  }

  return(genomic_distribution)

}

#' Plot Genomic Distribution
#'
#' Plot genomic distribution of TSSs or TSRs.
#'
#' @param genomic_distribution Dataframe of TSS or TSR genomic distributions from tsr_genomic_distribution
#'
#' @details
#' This plotting function will create a stacked barplot of genomic feature types for each sample.
#' The underlying data used for plotting comes from the 'genomic_distribution' function.
#' Genomic features include exons, introns, intergenic, downstream and promoter regions.
#' The promoter region is user defined when annotating the TSSs or TSRs.
#'
#' @return ggplot2 object with TSS or TSR genomic distribution plot
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
#' genomic_dist <- genomic_distribution(tsre_exp, data_type="tss")
#' plot_genomic_distribution(genomic_dist)
#'
#' @seealso \code{\link{annotate_features}} to annotate TSSs or TSRs to genomic features.
#'   \code{\link{genomic_distribution}} to prepare the annotated features for plotting.
#'
#' @rdname plot_genomic_distribution-function
#' @export 

plot_genomic_distribution <- function(
  genomic_distribution
) {

  ## Check inputs.
  assert_that(is(genomic_distribution, "DataFrame"))
  
  ## Pull out information from DataFrame.
  genomic_dist <- as_tibble(genomic_distribution, .name_repair="unique")

  ## Plot the genomic distribution.
  p <- ggplot(genomic_dist, aes(x=.data$samples, y=.data$count, fill=fct_rev(.data$simple_annotations))) +
    geom_col(position="fill") +
    coord_flip() +
    ylab("Fraction") +
    theme_bw() +
    theme(
      axis.title.y=element_blank(),
      panel.grid=element_blank()
    )

  if (metadata(genomic_distribution)$groupings) p <- p + facet_grid(fct_rev(factor(grouping)) ~ .)

  return(p)
}
