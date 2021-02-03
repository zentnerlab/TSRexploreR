#' TSR Metrics
#'
#' @description
#' Calculate various TSR Metrics including:
#'   IQR width, peak balance and shape index.
#'
#' @inheritParams common_params
#' @param iqr_lower Lower IQR cutoff value.
#' @param iqr_upper Upper IQR cutoff value.
#'
#' @description
#' TSR shape descriptors such as peaked versus broad can enrich for
#'   different functional classes of genes.
#' TSRexploreR allows for the calculation of various common TSR shape indicators
#'   such as Interquartile width (IQR), shape index (broad versus peaked), and peak balance.
#' 'iqr_lower' and 'iqr_upper' determine the upper and lower quantile bounds for calculation.
#'
#' @return TSRexplorer object with TSR metrics added to
#'   TSSs and TSRs.
#'
#' @examples
#' library("magrittr")
#' data(TSSs)
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer %>%
#'   format_counts(data_type="tss") %>%
#'   tss_clustering(threshold=3) %>%
#'   associate_with_tsr
#' tsr_metrics(tsre)
#'
#' @export

tsr_metrics <- function(
  experiment,
  iqr_lower=0.10,
  iqr_upper=0.90
) {

  ## Input Checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.numeric(iqr_lower) && (iqr_lower > 0 & iqr_lower < 1))
  assert_that(is.numeric(iqr_upper) && (iqr_upper > 0 & iqr_upper < 1))
  assert_that(iqr_lower < iqr_upper)

  ## Get samples from TSRexploreR object.
  select_samples <- experiment %>%
    extract_counts("tss", "all") %>%
    rbindlist(idcol="sample")

  keys <- c("sample", "TSR_FHASH")
  setkeyv(select_samples, keys)

  ## Calculate shape index.
  si <- shape_index(select_samples)
  setkeyv(si, keys)
  select_samples <- merge(select_samples, si, all.x=TRUE)

  ## Calculate peak concentration.
  #pc <- peak_concentration(select_samples)
  #setkeyv(pc, keys)
  #select_samples <- merge(select_samples, pc, all.x=TRUE)

  ## Calculate peak balance.
  pb <- peak_balance(select_samples)
  setkeyv(pb, keys)
  select_samples <- merge(select_samples, pb, all.x=TRUE)

  ## Calculate IQR.
  iqr <- iq_range(select_samples, iqr_upper, iqr_lower)
  setkeyv(iqr, keys)
  select_samples <- merge(select_samples, iqr, all.x=TRUE)

  ## Add metrics back to TSRs.
  tsr_names <- select_samples[["tsr_sample"]] %>%
    unique %>%
    discard(~ is.na(.))
  
  select_TSRs <- experiment %>%
    extract_counts("tsr", tsr_names) %>%
    rbindlist(idcol="tsr_sample")
  setnames(select_TSRs, old="FHASH", new="TSR_FHASH")

  tsr_metrics <- select_samples[
    !is.na(TSR_FHASH),
    .(tsr_sample, TSR_FHASH, shape_index, shape_class, peak_balance,
    iqr_min, iqr_max, iqr_width, iqr_coords)
  ]
  tsr_metrics <- unique(tsr_metrics)  

  tsr_keys <- c("tsr_sample", "TSR_FHASH")
  setkeyv(select_TSRs, tsr_keys)
  setkeyv(tsr_metrics, tsr_keys)

  select_TSRs <- merge(select_TSRs, tsr_metrics, all.x=TRUE)

  ## Add metrics back to TSRexploreR object.
  select_TSRs <- split(select_TSRs, select_TSRs$tsr_sample)
  walk(select_TSRs, function(x) {
    x[, tsr_sample := NULL]
    setnames(x, old="TSR_FHASH", new="FHASH")
    x <- as_granges(x)
    x <- as.data.table(x)
    return(x)
  })

  select_samples <- split(select_samples, select_samples$sample)
  walk(select_samples, function(x) {
    x[, sample := NULL]
    x <- as_granges(x)
    x <- as.data.table(x)
    return(x)
  })

  experiment <- set_count_slot(
    experiment, select_samples,
    "counts", "tss", "raw"
  )
  experiment <- set_count_slot(
    experiment, select_TSRs,
    "counts", "tsr", "raw"
  )
  return(experiment)
}

#' Shape Index
#'
#' Calculate shape index.
#'
#' @param tss_table data.table of TSSs

shape_index <- function(tss_table) {
  
  ## Calculate shape index.
  si_results <- tss_table[
    !is.na(TSR_FHASH) & score > 1,
    .(shape_index=((score / tsr_score) * log2(score / tsr_score))),
    by=.(sample, TSR_FHASH)
  ][,
    .(shape_index=2 + sum(shape_index)), by=.(sample, TSR_FHASH)
  ]

  ## Classify TSRs based on shape index.
  si_results[, shape_class := ifelse(shape_index < -1, "broad", "peaked")]
  
  ## Return results.
  return(copy(si_results)) 
}

#' Peak Concentration
#'
#' Calculate peak concentration.
#'
#' @param tss_table data.table of TSSs

peak_concentration <- function(tss_table) {
  
  ## Calculate peak concentration.
  pc_results <- tss_table[
    !is.na(TSR_FHASH) & score > 1,
    .(peak_concentration=log2((max(score) / sum(score)) * max(score))),
    by=.(sample, TSR_FHASH)
  ]

  ## Return peak concentration results.
  return(copy(pc_results))
}

#' Peak Balance
#'
#' Calculate peak balance.
#'
#' @param tss_table data.table of TSSs

peak_balance <- function(tss_table) {

  ## Get TSS position relative to TSR midpoints.
  tss_position <- tss_table[
    !is.na(TSR_FHASH) & score > 1,
    .(score, tsr_score, tss_pos=ifelse(
      strand == "+",
      start - median(range(start)),
      median(range(start)) - start
    )),
    by=.(sample, TSR_FHASH)
  ]

  ## Calculate peak balance.
  pb_results <- tss_position[,
    .(peak_balance=sum((score / sum(score)) * tss_pos)),
    by=.(sample, TSR_FHASH)
  ]

  ## Return peak balance values.
  return(copy(pb_results))

}

#' Interquartile Range
#'
#' Calculate IQR.
#'
#' @param tss_table data.table of TSSs
#' @inheritParams tsr_metrics

iq_range <- function(tss_table, iqr_upper, iqr_lower) {
  
  ## Get TSS positions relative to TSR midpoints.
  tss_position <- tss_table[
    !is.na(TSR_FHASH) & score > 1,
    .(score, seqnames, start, strand, tss_pos=ifelse(
      strand == "+",
      start - median(range(start)),
      median(range(start)) - start
    )),
    by=.(sample, TSR_FHASH)
  ]

  ## Calculate IQR.
  tss_position <- tss_position[order(sample, TSR_FHASH, tss_pos)]

  tss_position[,
    cum_sum := cumsum(score),
    by=.(sample, TSR_FHASH)
  ][,
    ecdf := ecdf_function(cum_sum),
    by=.(sample, TSR_FHASH)
  ]

  iqr_results <- tss_position[
    dplyr::between(ecdf, iqr_lower, iqr_upper)
  ][,
    .(iqr_min=min(cum_sum),
    iqr_max=max(cum_sum),
    iqr_width=max(tss_pos) - min(tss_pos),
    iqr_coords=str_c(seqnames, min(start), max(start), strand, sep=":")),
    by=.(sample, TSR_FHASH)
  ]
  iqr_results <- unique(iqr_results)

  ## Return IQR results.
  return(copy(iqr_results))

}

#' ECDF Function
#'
#' Return ECDF values.
#'
#' @param tss_vector Vector of values to compute ECDF values

ecdf_function <- function(tss_vector) {
  ecdf_func <- ecdf(tss_vector)
  ecdf_values <- ecdf_func(tss_vector)
  return(ecdf_values)
}
