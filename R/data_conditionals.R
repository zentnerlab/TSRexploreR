#' Data Conditionals
#'
#' @description
#' Apply strategies for filtering, ordering, quantiling, and/or grouping of data.
#'
#' @importFrom dplyr dense_rank ntile desc
#'
#' @param signal_data TSS or TSR data
#' @param quantile_by Continuous metric to calculate quantiles (qq??)
#' @param n_quantiles Number of quantiles to calculate for continuous metric
#' @param quantile_samples Samples to use when setting quantiles
#' @param quantile_group Group these features and quantile based on aggregate mean of feature
#' @param quantile_direction Order quantiles in 'descending' or 'ascending' order of values
#' @param order_by Metric to order data by
#' @param order_group Features wil be aggregated by this group before ordering
#' @param order_direction Whether the values should be ordered in
#'   'ascending' or 'descending' order
#' @param order_samples Names of samples to order by
#' @param filters Logical string to subset/filter data by
#' @param grouping If quantiles not set split data by categorical variable
#'
#' @details
#' TSSs and TSRs with different features, such as shape or sequence,
#'   have been connected with distinct biological functions.
#' Because of this, it may be desired to analyze certain subsets of TSSs or TSRs,
#'   or split the data based on various grouping strategies.
#' This function extends the flexibility of various other functions by adding the
#'   ability to filter, quantile, order, and/or group data prior to downstream analysis.
#'
#' 'filters' gives you the ability to filter TSSs or TSRs based on any column stored
#'   in the data.
#' The filters should be a character in standard R logical statements.
#' For example, if you wanted to keep only TSSs from chromosome I with a score greater than
#'   10, the filter would look like the following: "seqnames == 'chrI' & score > 10".
#'
#' The series of arguments for quantiling allows one to split the data based on the
#'   quantile of any column with numeric data.
#' 'quantile_by' is the column to use for calculating quantiles,
#'   and 'n_quantiles' specifies the number of quantiles to break the data into.
#' Whether the quantiles are calculated in ascending or descending order
#'   can be set with 'quantile_direction'.
#' If you want to use only a subset of samples to calculate the quantiles,
#'   they can be specified with 'quantile_samples'.
#' Finally, the data can first be split into subgroups before quantiling,
#'   such as TSR shape class.
#'
#' The ordering arguments can be used to control the order at which the TSSs or
#'   TSRs are plotted.
#' 'order_by' specifies the numeric column used to generate the order,
#'   and whether the order is ascending or descending is controlled by 'order_direction'.
#' 'order_group' can be used to first split the data based on some categorical value,
#'   such as TSR shape class.
#'
#' If quantiles are not specified, 'grouping' may be used to split the data based on
#'   a categorical value such as TSR shape class.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' tsre_exp <- tss_clustering(tsre_exp)
#' tsre_exp <- tsr_metrics(tsre_exp)
#' conditions <- list(order_by = "score", grouping = "shape_class")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "tsrexplorer")
#' seqs <- tss_sequences(
#'   tsre_exp, genome_assembly = assembly, threshold = 10,
#'   dominant = TRUE, data_conditions = conditions
#' )
#' plot_sequence_logo(seqs)
#'
#' @return data.table with the results from the various data conditionals applied.
#'
#' @rdname group_data-function
#' @export

group_data <- function(
  signal_data, filters = NA,
  order_by = NA, order_direction = "descending",
  order_samples = NA, order_group = NA,
  quantile_by = NA, n_quantiles = NA,
  quantile_samples = NA, quantile_group = NA,
  quantile_direction = "descending",
  grouping = NA
) {
  
  ## Filter the data if required  
  if (!is.na(filters)) {
    signal_data <- map(signal_data, function(x) {
      x <- subset(x, eval(parse(text = filters)))
      return(x)
    })
  }

  ## Order the data for plotting if required.
  if (!is.na(order_by)) {
    if (is.na(order_group)) order_group <- "FHASH"

    order_dataset <- rbindlist(signal_data, idcol = "sample")[,
      c("sample", ..order_by, ..order_group)
    ]

    order_dataset <- dcast(
      order_dataset, as.formula(str_c("sample ~ ", order_group)),
      value.var = order_by, fun.aggregate = mean,
      fill = ifelse(order_direction == "descending", -Inf, Inf)
    )

    order_dataset <- melt(
      order_dataset, id.vars = "sample",
      variable.name = order_group, value.name = order_by
    )

    if (!is.na(order_samples)) {
      order_dataset <- order_dataset[sample %in% order_samples]
    }

    order_dataset <- order_dataset[, .(order_by = mean(get(order_by))), by = order_group]
    setnames(order_dataset, old = "order_by", new = order_by)

    if (order_direction == "ascending") {
      order_dataset[, plot_order := dense_rank(desc(order_dataset[[order_by]]))]
    } else if (order_direction == "descending") {
      order_dataset[, plot_order := dense_rank(order_dataset[[order_by]])]
    }

    order_dataset[, c(order_by) := NULL]

    signal_data <- map(signal_data, function(x) {
        x <- merge(x, order_dataset, by = order_group, all.x = TRUE)
        return(x)
      })
  }

  ## Quantile the metric if required.
  if (!is.na(quantile_by)) {
    if (is.na(quantile_group)) quantile_group <- "FHASH"
  
    quantile_dataset <- rbindlist(signal_data, idcol = "sample")[,
      c("sample", ..quantile_by, ..quantile_group)
    ]

    quantile_dataset <- dcast(
      quantile_dataset, as.formula(str_c("sample ~ ", quantile_group)),
      value.var = quantile_by, fun.aggregate = mean,
      fill = ifelse(quantile_direction == "descending", -Inf, Inf)
    )

    quantile_dataset <- melt(
      quantile_dataset, id.vars = "sample",
      variable.name = quantile_group, value.name = quantile_by
    )

    if (!is.na(quantile_samples)) {
      quantile_dataset <- quantile_dataset[sample %in% quantile_samples]
    }

    quantile_dataset <- quantile_dataset[, .(quantile_by = mean(get(quantile_by))), by = quantile_group]
    setnames(quantile_dataset, old = "quantile_by", new = quantile_by)
    
    if (quantile_direction == "ascending") {
      quantile_dataset[, grouping := ntile(desc(quantile_dataset[[quantile_by]]), n_quantiles)]
    } else if (quantile_direction == "descending") {
      quantile_dataset[, grouping := ntile(quantile_dataset[[quantile_by]], n_quantiles)]
    }

    quantile_dataset[, c(quantile_by) := NULL]

    signal_data <- map(signal_data, function(x) {
      x <- merge(x, quantile_dataset, by = quantile_group, all.x = TRUE)
      return(x)
    })
  }

  ## If not using quantiles, split by a categorical value.
  if (!is.na(grouping) && is.na(quantile_by)) {
    signal_data <- map(signal_data, function(x) {
      x <- x[!is.na(x[[grouping]])]

      x[, grouping := x[[grouping]]]

      return(x)
    })
  }

  ## Return the signal data.
  return(signal_data)
}

#' Preliminary Filter
#'
#' Preliminary filtering of data
#'
#' @param signal_data TSS or TSR data
#' @param dominant Whether to retain dominant TSS/TSR
#' @param threshold Raw count threshold for a TSS or TSR to be considered
#'
#' @rdname preliminary_filter-function
#' @export

preliminary_filter <- function(signal_data, dominant, threshold) {
  
  ## Retain only dominant features if required.
  if (dominant) {
    signal_data <- map(signal_data, function(x) {
      x <- x[(dominant)]
      return(x)
    })
  }

  ## Apply a threshold to score if required.
  if (!is.na(threshold)) {
    signal_data <- map(signal_data, function(x) {
      x <- x[score >= threshold]
      return(x)
    })
  }

  ## Return signal data.
  return(signal_data)
}
