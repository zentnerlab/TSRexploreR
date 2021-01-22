#' Data Conditions Input
#'
#' @param data_filters Filter data by columns.
#' @param data_ordering Ordering object with order settings.
#'   See ?ordering for more information.
#' @param data_quantiling Quantiling object with quantile settings.
#'   See ?quantiling for more information.
#' @param data_grouping If quantiles not set,
#'   split data by the specified categorical variable.
#'
#' @export

conditionals <- function(
  data_filters=NULL,
  data_ordering=ordering(),
  data_quantiling=quantiling(),
  data_grouping=NULL
) {

  ## Check inputs.
  assert_that(
    is.null(data_ordering) || (
      is.list(data_ordering) && 
      has_attr(data_ordering, "condition_type") &&
      attributes(data_ordering)$condition_type == "ordering"
    )
  )
  assert_that(
    is.null(data_quantiling) || (
      is.list(data_quantiling) &&
      has_attr(data_quantiling, "condition_type") &&
      attributes(data_quantiling)$condition_type == "quantiling"
    )
  )

  ## Return a quosure of the filter.
  data_filters <- enquo(data_filters)
  if (quo_is_null(data_filters)) data_filters <- NULL

  ## For data grouping return either NULL or a character.
  if (quo_is_null(enquo(data_grouping))) {
    data_grouping <- NULL
  } else {
    data_grouping <- as.character(ensym(data_grouping))
  }

  ## Return a list of conditions.
  conds <- list(
    filters=data_filters,
    ordering=data_ordering,
    quantiling=data_quantiling,
    grouping=data_grouping
  )
  attr(conds, "data_conditions") <- "data_conditions"

  return(conds)  

}

#' Ordering
#'
#' @param ... Variables to order by.
#'   Wrap varaible name in desc() for descending order (like in dplyr::arrange).
#' @param samples Names of samples to order by aggregate score.
#' @param If more than one sample is selected feature values
#'   are aggregated using this function.
#'
#' @export

ordering <- function(
  ...,
  .samples=NULL,
  .aggr_fun=mean
) {

  ord <- enquos(...)

  if (length(ord) == 0) {
    ord <- NULL
  } else {
    ord <- list(
      ordering=ord,
      samples=.samples,
      aggr_fun=.aggr_fun
    )
    attr(ord, "condition_type") <- "ordering"
  }

  return(ord)
}

#' Quantiling
#'
#' @param by Continuous metric for calculating quantiles.
#' @param n Number of quantiles to calculate for continuous metric.
#' @param samples Samples to use when setting quantiles.
#' @param descending Order quantiles in descending order.
#' @param aggr_fun Function to aggregate scores if more than one sample selected.
#'
#' @export

quantiling <- function(
    by=NULL,
    n=NULL,
    samples=NULL,
    descending=TRUE,
    aggr_fun=mean
){

  ## Check inputs.
  assert_that(is.null(n) || is.count(n))
  assert_that(is.null(samples) || is.character(samples))
  assert_that(is.flag(descending))
  assert_that(is.function(aggr_fun))

  ## Quosure for by column.
  by <- enquo(by)

  ## Return NULL or list of values.
  if (quo_is_null(by)) {
    quantiling <- NULL
  } else {
    quantiling <- list(
      by=quo_text(by),
      n=n,
      samples=samples,
      descending=descending,
      aggr_fun=aggr_fun
    )
    attr(quantiling, "condition_type") <- "quantiling"
  }

  return(quantiling)

}

#' Data Conditionals
#'
#' @description
#' Apply strategies for filtering, ordering, quantiling, and/or grouping of data.
#'
#' @importFrom dplyr dense_rank ntile desc
#'
#' @param signal_data TSS or TSR data.
#' @param data_conditions List of conditions.
#'
#' @details
#' It may be desirable to analyze certain subsets of TSSs or TSRs, or split the data 
#' based on various categorical variables. This function extends the flexibility of 
#' various other functions by adding the ability to filter, quantile, order, and/or 
#' group data prior to downstream analysis.
#'
#' 'filters' gives you the ability to filter TSSs or TSRs based on any column stored
#' in the data. The filters should be a character in standard R logical statements. 
#' For example, if you wanted to keep only TSSs from chromosome I with a score greater than
#' 10, the filter would be: "seqnames == 'chrI' & score > 10".
#'
#' The series of arguments for quantiling allows one to split the data based on the
#'   quantile of any column with numeric data. 'quantile_by' is the column to use 
#'   for calculating quantiles, and 'n_quantiles' specifies the number of quantiles 
#'   to break the data into. Whether the quantiles are calculated in ascending or 
#'   descending order can be set with 'quantile_direction'. If you want to use only 
#'   a subset of samples to calculate the quantiles, they can be specified with 'quantile_samples'. 
#'   Finally, the data can first be split into subgroups before quantiling,
#'   such as TSR shape class.
#'
#' The ordering arguments can be used to control the order in which the TSSs or
#' TSRs are plotted. 'order_by' specifies the numeric column used to generate the order,
#' and whether the order is ascending or descending is controlled by 'order_direction'.
#' 'order_group' can be used to first split the data based on some categorical value,
#' such as TSR shape class.
#'
#' If quantiles are not specified, 'grouping' may be used to split the data based on
#' a categorical value such as TSR shape class.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' exp <- tsr_explorer(TSSs)
#' exp <- format_counts(exp, data_type="tss")
#' exp <- tss_clustering(exp)
#' exp <- tsr_metrics(exp)
#' conditionals <- list(order_by="score", grouping="shape_class")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#' seqs <- tss_sequences(
#'   exp, genome_assembly=assembly, threshold=10,
#'   dominant=TRUE, data_conditions=conditionals
#' )
#' plot_sequence_logo(seqs)
#'
#' @return data.table with the results from the various data conditionals applied.
#'
#' @rdname group_data-function
#' @export

condition_data <- function(
  signal_data,
  data_conditions
) {

  ## Input checks.
  assert_that(
    is.null(data_conditions) || (
      is.list(data_conditions) &&
      has_attr(data_conditions, "data_conditions")
    )
  )

  ## Filter the data if required.
  if (!is.null(data_conditions$filters)) {
    walk(signal_data, setDT)
    signal_data <- .filter_data(signal_data, data_conditions$filters)
  }

  ## Ordering the data if required.
  if (!is.null(data_conditions$ordering)) {
    walk(signal_data, setDT)
    signal_data <- .order_data(signal_data, data_conditions$ordering)
  }

  ## Quantiling the data if required.
  if (!is.null(data_conditions$quantiling)) {
    walk(signal_data, setDT)
    signal_data <- .quantile_data(signal_data, data_conditions$quantiling)
  }

  ## Grouping the data if required.
  ## Modifies the tables in place.
  if (!is.null(data_conditions$grouping)) {
    walk(signal_data, setDT)
    .group_data(signal_data, data_conditions$grouping)
  }

  ## Return the signal data.
  return(signal_data)
}

#' Filtering
#'
#' @param signal_data List of TSSs or TSRs
#' @param filters Quosure containing filters

.filter_data <- function(
  signal_data,
  filters
) {
  signal_data <- map(signal_data, ~dplyr::filter(.x, !!filters))
  return(signal_data)
}

#' Ordering
#'
#' @param signal_data List of TSSs or TSRs
#' @param ordering List of ordering parameters

.order_data <- function(
  signal_data,
  ordering
) {

  ## Create consensus ranges with aggregated scores.
  signal_data <- rbindlist(signal_data, idcol="sample")
  setkey(signal_data, seqnames, strand, start, end)

  ord <- .create_consensus(signal_data, ordering)
  
  ## Arrange the aggregated features and add numeric ordering.
  ord <- dplyr::arrange(ord, !!!ordering$ordering)
  setDT(ord)
  ord[, row_order := .I]
  ord <- ord[, .(CFHASH, row_order)]
  ord[, c("seqnames", "start", "end", "strand") :=
    tstrsplit(CFHASH, split=":", fixed=TRUE, type.convert=TRUE)
  ]
  ord[, CFHASH := NULL]

  ## Merge the row order into the original data.
  setkey(ord, seqnames, strand, start, end)
  signal_data <- foverlaps(ord, signal_data)
  signal_data[, c("i.start", "i.end") := NULL]

  ## Return the signal data.
  signal_data <- split(signal_data, by="sample", keep.by=FALSE)
  return(signal_data)
  
}

#' Quantiling
#'
#' @param signal_data list of TSSs or TSRs
#' @param quantiling list of quantiling parameters

.quantile_data <- function(
  signal_data,
  quantiling
) {

  ## Create consensus ranges with aggregated scores.
  signal_data <- rbindlist(signal_data, idcol="sample")
  setkey(signal_data, seqnames, strand, start, end)

  ord <- .create_consensus(signal_data, quantiling)

  ## quantile the aggregated features.
  if (quantiling$descending) {
    ord[, row_quantile := ntile(desc(ord[[quantiling$by]]), n=quantiling$n)]
  } else {
    ord[, row_quantile := ntile(ord[[quantiling$by]], n=quantiling$n)]
  }

  ord <- ord[, .(CFHASH, row_quantile)]
  ord[, c("seqnames", "start", "end", "strand") :=
    tstrsplit(CFHASH, split=":", fixed=TRUE, type.convert=TRUE)
  ]
  ord[, CFHASH := NULL]

  ## Merge the row order into the original data.
  setkey(ord, seqnames, strand, start, end)
  signal_data <- foverlaps(ord, signal_data)
  signal_data[, c("i.start", "i.end") := NULL]

  ## Return the signal data.
  signal_data <- split(signal_data, by="sample", keep.by=FALSE)
  return(signal_data)

}

#' Group Data
#'
#' @param signal_data TSS or TSR data
#' @param grouping Grouping variable

.group_data <- function(
  signal_data,
  grouping
) {
  walk(signal_data, function(x) {
    x[, row_groups := x[[grouping]]]
  })
}

#' Create Consensus Ranges
#'
#' @param signal_data TSS or TSR data
#' @param conditionals Either ordering or quantiling conditions

.create_consensus <- function(
  signal_data,
  conditionals
) {

  ## Reduce ranges of samples into consensus ranges.
  cranges <- signal_data %>%
    as_granges %>%
    reduce_ranges_directed %>%
    as.data.table(key=c("seqnames", "strand", "start", "end"))
  cranges[, CFHASH := str_c(seqnames, start, end, strand, sep=":")]

  ## Filter out samples not being used to calculate ordering.
  if (!is.null(conditionals$samples)) {
    ord <- signal_data[samples %in% conditionals$samples]
  } else {
    ord <- signal_data
  }
  setkey(ord, seqnames, strand, start, end)

  ## Merge consensus ranges into ranges used for ordering.
  ord <- foverlaps(cranges, ord)
  ord[, c(
    "i.start", "i.end", "i.width", "seqnames",
    "strand", "start", "end", "FHASH"
  ) := NULL]

  ## Aggregate numeric features.
  numeric_cols <- colnames(ord)[sapply(ord, is.numeric)]

  ord <- ord[,
    lapply(.SD, conditionals$aggr_fun),
    by=CFHASH,
    .SDcol=numeric_cols
  ]

  return(ord)
}

#' Preliminary Filter
#'
#' Preliminary filtering of data
#'
#' @param signal_data TSS or TSR data.
#' @param dominant Whether to use dominant TSS/TSR.
#' @param threshold Raw count threshold for a TSS or TSR to be considered.
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
  if (!is.null(threshold)) {
    signal_data <- map(signal_data, function(x) {
      x <- x[score >= threshold]
      return(x)
    })
  }

  ## Return signal data.
  return(signal_data)
}
