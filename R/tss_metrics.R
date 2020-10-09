
#' Mark Dominant
#'
#' @description
#' Mark TSSs as dominant TSS per TSR or gene, or TSRs as dominant per gene. 
#'
#' @param experiment tsrexplorer object with annotated TSSs/TSRs
#' @param data_type Either 'tss' or 'tsr'
#' @param threshold Read threshold for TSS/TSRs
#' @param mark_per By default marks dominant TSR per gene, and dominant TSS per TSR.
#'   TSSs can also be set per as dominant TSS per 'gene'.
#'
#' @details
#' This function marks which TSSs are dominant per TSR or gene,
#'   or which TSR is dominant per gene.
#' Analysis of dominant features may help to cut through the noise to get
#'   information such as the primary 5' UTR, sequence features associated with the
#'   the strongest TSS, and other related questions.
#'
#' Setting a 'threshold' will only mark a TSS or TSR as dominant if their score
#'    is greater than or equal to the threshold.
#'
#' 'mark_per' controls the behavior of the function.
#' For TSSs 'default' will mark dominant TSS per TSR, and for TSRs the dominant
#'   TSR per gene is marked.
#' for TSSs, 'gene' can also be specified, which will mark the dominant TSS per gene.  
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#' tsre_exp <- associate_with_tsr(tsre_exp, sample_list = list(
#'   "S288C_WT_1" = "S288C_WT_1", "S288C_WT_2" = "S288C_WT_2", "S288C_WT_3" = "S288C_WT_3",
#'   "S288C_D_1" = "S288C_D_1", "S288C_D_2" = "S288C_D_2", "S288C_D_3" = "S288C_D_3"
#' ))
#' tsre_exp <- mark_dominant(tsre_exp, data_type = "tss")
#'
#' @return tsr exlorer object with dominant status added to TSSs or TSRs.
#'
#' @seealso
#' \code{\link{associate_wth_tsr}} to associate TSSs with TSRs prior to marking
#'   dominant TSS per TSR.
#'
#' @rdname mark_dominant-function
#' @export

mark_dominant <- function(
  experiment,
  data_type = c("tss", "tsr"),
  threshold = 1,
  mark_per = "default"
) {

  ## Check inputs.
  if (!is(experiment, "tsr_explorer")) stop("'experiment' must be a tsr explorer object")
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr"))
  if (
    !is.na(threshold) && !is(threshold, "numeric") ||
    threshold %% 1 != 0 || threshold < 1
  ) {
    stop("threshold must be a positive integer greater than or equal to 1")
  }
  mark_per <- match.arg(str_to_lower(mark_per), c("default", "gene"))
  
  ## Select samples.
  select_samples <- extract_counts(experiment, data_type, "all")

  ## Mark dominant TSS/TSR per gene if requested.
  if (data_type == "tsr" | (data_type == "tss" & mark_per == "gene")) {
    dominant <- map(select_samples, function(x) {
      x[,
        dominant := (
          score == max(score) &
          !simple_annotations %in% c("Downstream", "Intergenic") &
          score >= threshold
        ),
        by = eval(ifelse(
          experiment@settings$annotation[, feature_type] == "transcript",
          "transcriptId", "geneId"
        ))
      ]

      return(x)
    })
  
  ## Mark the dominant TSS per TSR if requested.
  } else if (data_type == "tss" & mark_per == "default") {
    dominant <- map(select_samples, function(x) {
      x[,
        dominant := (
          !is.na(score) &
          score == max(score) &
          score >= threshold
        ),
        by = TSR_FID
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

#' Max UTR Length
#'
#' Get TSS with furthest distance
#'
#' @param experiment tsrexplorer object with annotated TSSs
#' @param samples Either 'all' or names of sample to analyze
#' @param threshold Number of reads required for each TSS
#' @param max_upstream Max upstream distance of TSS to consider
#' @param max_downstream Max downstream distance of TSS to consider
#' @param feature_type Feature type used when finding distance to TSS ("geneId", "transcriptId")
#' @param quantiles Number of quantiles to break data into.
#'
#' @return tibble with max UTR length for features
#'
#' @rdname max_utr-function
#'
#' @export

max_utr <- function(
  experiment,
  samples = "all",
  threshold = 1,
  max_upstream = 1000,
  max_downstream = 100,
  feature_type = c("geneId", "transcriptId"),
  quantiles = NA
) {
  ## Grab selected samples.
  max_utr <- experiment %>%
    extract_counts("tss", samples) %>%
    bind_rows(.id = "sample") %>%
    as.data.table

  setnames(max_utr, old = feature_type, new = "feature")

  ## Filter data.
  max_utr <- max_utr[
    score >= threshold &
    dplyr::between(distanceToTSS, -max_upstream, max_downstream),
    .(sample, distanceToTSS, feature, score)
  ]

  ## Get TSS with minimum distance to start codon.
  max_utr <- max_utr[,
    .SD[which.min(distanceToTSS)],
    by = .(sample, feature)
  ]

  ## Add quantiles if requested.
  if (!is.na(quantiles)) {
    max_utr[, ntile := ntile(score, quantiles), by = sample]
  }

  ## Return DataFrame.
  max_utr <- DataFrame(max_utr)
  metadata(max_utr)$quantiles <- quantiles
  metadata(max_utr)$max_upstream <- max_upstream
  metadata(max_utr)$max_downstream <- max_downstream
  metadata(max_utr)$threshold <- threshold

  return(max_utr)
}

#' Plot Max UTR Length
#'
#' Plot TSS with furthest distance
#'
#' @param max_utr tibble of max UTRs output by max_utr
#' @param upstream Bases upstream to extend average to
#' @param downstream Bases downstream to extend average to
#' @param ncol Number of columns to plot the data to
#' @param consider_score Whether the score of the TSS should be
#' considered when plotting
#' @param ... Arguments passed to geom_density
#'
#' @return ggplot2 object of max UTR length average
#'
#' @rdname plot_max_utr-function
#'
#' @export

plot_max_utr <- function(
  max_utr, upstream = 1000, downstream = 100,
  ncol = 1, consider_score = FALSE, ...
) {

  ## Grab some info from DataFrame.
  quantiles <- metadata(max_utr)$quantiles

  ## Prepare data for plotting.
  utr_plot <- as.data.table(max_utr)
  
  if (!is.na(quantiles)) {
    utr_plot[, ntile := fct_rev(factor(ntile))]
  }

  ## Format data if score should be considered.
  if (consider_score) {
    utr_plot <- utr_plot[rep(seq_length(.N), score)]
  }

  ## Plot max UTR length detected.
  p <- ggplot(utr_plot, aes(x = .data$distanceToTSS)) +
    geom_density(fill = "#431352", color = "#431352")+#, ...) +
    xlim(-upstream, downstream) +
    theme_bw() +
    labs(
      x = "Max UTR Length",
      y = "Density"
    ) +
    geom_vline(xintercept = 0, lty = 2)

  if (!is.na(quantiles)) {
    p <- p + facet_grid(ntile ~ sample)
  } else {
    p <- p + facet_wrap(~ sample, ncol = ncol)
  }

  return(p)
}
