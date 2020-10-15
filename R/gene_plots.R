
#' Genes/Transcripts Detected
#'
#' @description
#' Get the number of genes or transcripts with an associated unique TSS or TSR
#'
#' @param experiment tsrexplorer object with annotated TSSs or TSRs
#' @param samples Either 'all' or a vector of sample names
#' @param data_type Whether TSSs or TSRs should be analyzed
#' @param threshold The number of raw reads required for a TSS or TSR to be considered
#' @param dominant Whether to consider only the dominant TSS or TSR
#' @param condition_data Apply conditions to data (supports filtering and quantiles/grouping)
#'
#' @details
#' This function will return either the number of genes or transcripts with
#'  an associated unique TSS or TSR.
#' Information on whether the feature has a promoter-proximal TSS or TSR is included
#'   in the output for plotting purposes.
#'
#' A set of functions to control data structure for plotting are included.
#' 'threshold' will define the minimum number of reads a TSS or TSR
#'  must have to be considered.
#' 'dominant' specifies whether only the dominant TSS or TSR is considered 
#'   from the 'mark_dominant' function.
#' For TSSs this can be either dominant per TSR or gene, and for TSRs
#'   it is just the dominant TSR per gene.
#' 'data_conditions' allows for the advanced filtering, ordering, and grouping
#'   of data.
#'
#' @return DataFrame of detected feature numbers.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data = annotation,
#'   data_type = "tss", feature_type = "transcript"
#' )
#' detected <- detect_features(tsre_exp, data_type = "tss")
#'
#' @seealso
#' \code{\link{annotate_features}} to annotate the TSSs and TSRs.
#' \code{\link{plot_detected_features}} to plot detected features.
#'
#' @rdname detect_features-function
#' @export

detect_features <- function(
  experiment,
  samples = "all",
  data_type = c("tss", "tsr", "tss_features", "tsr_features"),
  threshold = NA,
  dominant = FALSE,
  condition_data = NA
) {

  ## Check inputs.
  if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")
  if (!is(samples, "character")) stop("samples must be a character vecotor")
  data_type <- match.arg(str_to_lower(data_type), c("tss", "tsr", "tss_features", "tsr_features"))
  if (
    !is.na(threshold) && (!is(threshold, "numeric") ||
    threshold %% 1 != 0 || threshold < 1)
  ) {
    stop("threshold must be a positive integer greater than or equal to 1")
  }
  if (!is(dominant, "logical")) stop("dominant must be logical")
  if (all(!is.na(condition_data)) && !is(condition_data, "list")) {
    stop("condition_data must be a list of values")
  }
  
  ## Get sample data.
  sample_data <- extract_counts(experiment, data_type, samples)

  ## Initial sample processing.
  if (data_type %in% c("tss", "tsr") && (dominant | !is.na(threshold))) {
    sample_data <- preliminary_filter(sample_data, dominant, threshold)
  }

  if (data_type %in% c("tss_features", "tsr_features")) {
    sample_data <- map(sample_data, ~ .x[score >= threshold])
  }
  
  ## Apply data conditioning if set.
  if (data_type %in% c("tss", "tsr") && all(!is.na(condition_data))) {
    sample_data <- do.call(group_data, c(list(signal_data = sample_data), condition_data))
  }

  ## Rename feature column.
  walk(sample_data, function(x) {
    setnames(
      x, old = ifelse(
        experiment@settings$annotation[, feature_type] == "transcript",
        "transcriptId", "geneId"
      ),
      new = "feature"
    )
  })

  ## Get feature counts.
  groupings <- data_type %in% c("tss", "tsr") &&
    any(names(condition_data) %in% c("quantile_by", "grouping"))
  sample_data <- rbindlist(sample_data, idcol = "sample")

  if (data_type %in% c("tss", "tsr")) {
    if (groupings) {
      sample_data <- sample_data[,
        .(grouping, promoter = any(simple_annotations == "Promoter")),
        by = .(sample, feature)
      ][,
        .(with_promoter = sum(promoter), without_promoter = .N - sum(promoter), total = .N),
        by = .(sample, grouping)
      ]
    } else {
      sample_data <- sample_data[,
        .(promoter = any(simple_annotations == "Promoter")),
        by = .(sample, feature)
      ][,
        .(with_promoter = sum(promoter), without_promoter = .N - sum(promoter), total = .N),
        by = sample
      ]
    }
  } else if (data_type %in% c("tss_features", "tsr_features")) {
    sample_data <- sample_data[, .(count = uniqueN(feature)), by = sample]
  }

  ## Order samples if required.
  if (!all(samples == "all")) {
    sample_data[, samples := factor(samples, levels = samples)]
  }

  ## Create DataFrame to export.
  detected_features <- DataFrame(sample_data)

  metadata(detected_features)$threshold <- threshold
  metadata(detected_features)$data_type <- data_type
  metadata(detected_features)$feature_type <- experiment@settings$annotation[, feature_type]
  if (data_type %in% c("tss", "tsr")) {
    metadata(detected_features)$dominant <- dominant
    metadata(detected_features)$groupings <- groupings
  }

  return(detected_features)
}

#' Plot Detected Features
#'
#' @description
#' Plot number of features detected per sample.
#'
#' @importFrom stringr str_to_title
#'
#' @param detected_features DataFrame of detected feature counts from detect_features
#' @param ... Arguments passed to geom_col
#'
#' @details
#' This plotting function returns a stacked bar-plot showing the number of
#'   features detected with and without a promoter proximal TSS or TSR.
#' The information is first prepared using the 'detect_features' function.
#'
#' @return ggplot2 object of detected feature counts
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data = annotation,
#'   data_type = "tss", feature_type = "transcript"
#' )
#' detected <- detect_features(tsre_exp, data_type = "tss")
#' plot_detected_features(detected)
#'
#' @seealso
#' \code{\link{annotate_features}} to annotate the TSSs or TSRs.
#' \code{\link{detect_features}} to first detect feature numbers.
#'
#' @rdname plot_detected_features-function
#' @export

plot_detected_features <- function(
  detected_features,
  ...
) {

  ## Check inputs.
  if (!is(detected_features, "DataFrame")) stop("detected_features must be a DataFrame")
  
  ## Get some info from DataFrame.
  data_type <- metadata(detected_features)$data_type
  feature_type <- str_to_title(metadata(detected_features)$feature_type)
  if (data_type %in% c("tss", "tsr")) {
    grouping_status <- metadata(detected_features)$groupings
  }

  ## Prepare data for plotting.
  plot_data <- as.data.table(detected_features)

  if (data_type %in% c("tss", "tsr")) {
    plot_data[, total := NULL]  
    plot_data <- melt(plot_data,
      measure.vars = c("with_promoter", "without_promoter"),
      variable.name = "count_type",
      value.name = "feature_count"
    )
    plot_data[, count_type := fct_rev(factor(count_type))]
  }

  ## Plot data.
  if (data_type %in% c("tss", "tsr")) {
    p <- ggplot(plot_data, aes(x = .data$sample, y = .data$feature_count, fill = .data$count_type)) +
      geom_col(position = "stack", ...) +
      theme_bw() +
      ylim(c(0, NA)) +
      ylab(str_c(feature_type, "Count", sep = " ")) +
      xlab("Sample") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    if (grouping_status) {
      p <- p + facet_grid(fct_rev(factor(grouping)) ~ .)
    }
  } else if (data_type %in% c("tss_features", "tsr_features")) {
    p <- ggplot(plot_data, aes(x = .data$sample, y = .data$count, fill = .data$sample)) +
      geom_col(...) +
      theme_bw() +
      scale_fill_viridis_d() +
      ylim(c(0, NA)) +
      ylab(str_c(feature_type, "Count", sep = " ")) +
      xlab("Sample") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  }
    
  return(p)
}
