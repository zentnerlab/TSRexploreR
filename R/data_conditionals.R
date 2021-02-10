#' Data Conditions Input
#'
#' @description
#' This function allows  filtering, ordering, quantiling,
#'   and grouping of data for plotting.
#'
#' @param data_filters Logical statements by which to filter data.
#' @param data_ordering Order object with order settings.
#'   See ?ordering for more information.
#' @param data_quantiling Quantile object with quantile settings.
#'   See ?quantiling for more information.
#' @param data_grouping If quantiles not set,
#'   split data by the specified categorical variable.
#'
#' @details
#' It may be desirable to analyze certain subsets of TSSs or TSRs, or split the data 
#' based on various categorical variables. This function extends the flexibility of 
#' various other functions by adding the ability to filter, quantile, order, and/or 
#' group data prior to downstream analysis.
#'
#' 'data_filters' takes logical statements to filter TSSs or TSRs
#'   by any column stored in the data.
#' 'data_ordering' takes an 'ordering' object as input,
#'   which allows ordering of data by one or more columns.
#' 'data_quantiling' takes a 'quantiling' object as input,
#'   and will split plots by the given number of quantiles.
#' 'data_grouping' will split a plot by the given column if
#'   quantiling is not set.
#'
#' @return 'data_conditions' object for input to the corresponding
#'   'data_conditions' argument in select functions.
#'
#' @examples
#' data(TSSs)
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#'
#' exp <- TSSs[1] %>%
#'   tsr_explorer(genome_assembly=assembly) %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3) %>%
#'   associate_with_tsr %>%
#'   tsr_metrics
#'
#' # Sequence logo of TSSs from peaked TSRs
#' conditions <- conditionals(shape_class == "peaked")
#' plot_sequence_logo(exp, data_conditions=conditions)
#'
#' # Sequence color map sorted by descending TSS score
#' conditions <- conditionals(data_ordering=ordering(desc(score)))
#' plot_sequence_colormap(exp, data_conditions=conditions)
#'
#' # Sequence logos split by TSS score quantile
#' conditions <- conditionals(data_quantiling=quantiling(score, n=5))
#' plot_sequence_logo(exp, data_conditions=conditions)
#'
#' # Sequence logo split by TSR shape class
#' conditions <- conditionals(data_grouping=shape_class)
#' plot_sequence_logo(exp, data_conditions=conditions)
#'
#' @seealso
#' \code{\link{ordering}} for ordering information.
#' \code{\link{quantiling}} for quantiling information.
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
#' @description
#' This is a companion function to 'conditionals' that allows the
#'   specification of variables to order data by for plotting.
#'
#' @param ... Variables to order by.
#'   Wrap varaible name in desc() for descending order (like in dplyr::arrange).
#' @param .samples Names of samples to order by aggregate score.
#' @param .aggr_fun If more than one sample is selected feature values
#'   are aggregated using this function.
#'
#' @description
#' Columns to order by should be specified by names/symbols and not characters.
#' To sort by a variable in descending order, wrap the variable name in the 'desc' function.
#' By default ordering is calculated based on the aggregate value (based on '.aggr_fun')
#'   of the variable accross all samples.
#' If '.samples' is specified, only the given samples will be used for aggregate
#'   calculation and ordering.
#'
#' @return A list of ordering parameters to be passed to the 'data_ordering' argument
#'   of 'conditionals'.
#'
#' @examples
#' \donttest{
#' data(TSSs)
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#'
#' exp <- TSSs[1] %>%
#'   tsr_explorer(genome_assembly=assembly) %>%
#'   format_counts(data_type="tss")
#'
#' # Sequence color map sorted by descending score.
#' conditions <- conditionals(data_ordering=ordering(desc(score)))
#' plot_sequence_colormap(exp, data_conditions=conditions)
#' }
#'
#' @seealso
#' \code{\link{conditionals}} For more information on advanced data conditions.
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
#' @description
#' This is a companion function to 'conditionals' that allows the
#'   specification of a variable by which to split the data into 
#'   quantiles before plotting.
#'
#' @param by Continuous metric for calculating quantiles.
#' @param n Number of quantiles to calculate.
#' @param samples Samples to use when setting quantiles.
#' @param descending Order quantiles in descending order.
#' @param aggr_fun Function by which to aggregate scores if more than one sample selected.
#'
#' @details
#' Column to order by should be specified by names/symbol and not character.
#' By default ,quantiling is calculated based on the aggregate value (based on 'aggr_fun')
#'   of the variable across all samples.
#' If 'samples' is specified, only the given samples will be used for aggregate
#'   calculation and ordering.
#' 'descending' controls whether quantiling is calculated in descending (TRUE) or
#'   ascending (FALSE) order, and 'n' allows specification of number of quantiles.
#' 
#' @return A list of quantiling parameters to be passed to the 'data_quantiling' argument
#'   of 'conditionals'.
#'
#' @examples
#' \donttest{
#' data(TSSs)
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#'
#' exp <- TSSs[1] %>%
#'   tsr_explorer(genome_assembly=assembly) %>%
#'   format_counts(data_type="tss")
#'
#' # Sequence base color plot quantiled by score.
#' conditions <- conditionals(data_quantiling=quantiling(score, n=5))
#' plot_sequence_colormap(exp, data_conditions=conditions)
#' }
#'
#' @seealso
#' \code{\link{conditionals}} For more information on advanced data conditions.
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
#' @param signal_data List of TSSs or TSRs
#' @param quantiling List of quantiling parameters

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
