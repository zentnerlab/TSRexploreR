
#' TSR Metrics
#'
#' Calculate TSR Metrics.
#'
#' @param experiment tsrexplorer object
#'
#' @rdname tsr_metrics-function
#' @export

tsr_metrics <- function(experiment) {

  ## Input Checks.
  if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")

  ## Get samples from tsrexplorer object.
  select_samples <- experiment %>%
    extract_counts("tss", "all") %>%
    rbindlist(idcol = "sample")

  keys <- c("sample", "TSR_FID")
  setkeyv(select_samples, keys)

  ## Calculate shape index.
  si <- shape_index(select_samples)
  setkeyv(si, keys)
  select_samples <- merge(select_samples, si, all.x = TRUE)

  ## Calculate peak concentration.
  #pc <- peak_concentration(select_samples)
  #setkeyv(pc, keys)
  #select_samples <- merge(select_samples, pc, all.x = TRUE)

  ## Calculate peak balance.
  pb <- peak_balance(select_samples)
  setkeyv(pb, keys)
  select_samples <- merge(select_samples, pb, all.x = TRUE)

  ## Calculate IQR.
  iqr <- iq_range(select_samples)
  setkeyv(iqr, keys)
  select_samples <- merge(select_samples, iqr, all.x = TRUE)

  ## Add metrics back to TSRs.
  tsr_names <- select_samples[["tsr_sample"]] %>%
    unique %>%
    discard(~ is.na(.))
  
  select_TSRs <- experiment %>%
    extract_counts("tsr", tsr_names) %>%
    rbindlist(idcol = "tsr_sample")
  setnames(select_TSRs, old = c("FID", "FHASH"), new = c("TSR_FID", "TSR_FHASH"))

  tsr_metrics <- select_samples[
    !is.na(TSR_FID),
    .(tsr_sample, TSR_FID, TSR_FHASH, shape_index, shape_class, peak_balance,
    iqr_min, iqr_max, iqr_width, iqr_coords)
  ]
  tsr_metrics <- unique(tsr_metrics)  

  tsr_keys <- c("tsr_sample", "TSR_FID", "TSR_FHASH")
  setkeyv(select_TSRs, tsr_keys)
  setkeyv(tsr_metrics, tsr_keys)

  select_TSRs <- merge(select_TSRs, tsr_metrics, all.x = TRUE)

  ## Add metrics back to tsrexplorer object.
  select_TSRs <- split(select_TSRs, select_TSRs$tsr_sample)
  walk(select_TSRs, function(x) {
    x[, tsr_sample := NULL]
    setnames(x, old = c("TSR_FID", "TSR_FHASH"), new = c("FID", "FHASH"))
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
#'
#' @rdname shape_index-function
#' @export

shape_index <- function(tss_table) {
  
  ## Calculate shape index.
  si_results <- tss_table[
    !is.na(TSR_FID) & score > 1,
    .(shape_index = ((score / tsr_score) * log2(score / tsr_score))),
    by = .(sample, TSR_FID)
  ][,
    .(shape_index = 2 + sum(shape_index)), by = .(sample, TSR_FID)
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
#'
#' @rdname peak_concentration-function
#' @export

peak_concentration <- function(tss_table) {
  
  ## Calculate peak concentration.
  pc_results <- tss_table[
    !is.na(TSR_FID) & score > 1,
    .(peak_concentration = log2((max(score) / sum(score)) * max(score))),
    by = .(sample, TSR_FID)
  ]

  ## Return peak concentration results.
  return(copy(pc_results))
}

#' Peak Balance
#'
#' Calculate peak balance.
#'
#' @param tss_table data.table of TSSs
#'
#' @rdname peak_balance-function
#' @export

peak_balance <- function(tss_table) {

  ## Get TSS position relative to TSR midpoints.
  tss_position <- tss_table[
    !is.na(TSR_FID) & score > 1,
    .(score, tsr_score, tss_pos = ifelse(
      strand == "+",
      start - median(range(start)),
      median(range(start)) - start
    )),
    by = .(sample, TSR_FID)
  ]

  ## Calculate peak balance.
  pb_results <- tss_position[,
    .(peak_balance = sum((score / sum(score)) * tss_pos)),
    by = .(sample, TSR_FID)
  ]

  ## Return peak balance values.
  return(copy(pb_results))

}

#' Interquartile Range
#'
#' Calculate IQR.
#'
#' @param tss_table data.table of TSSs
#'
#' @rdname iq_range-function
#' @export

iq_range <- function(tss_table) {
  
  ## Get TSS positions relative to TSR midpoints.
  tss_position <- tss_table[
    !is.na(TSR_FID) & score > 1,
    .(score, seqnames, start, strand, tss_pos = ifelse(
      strand == "+",
      start - median(range(start)),
      median(range(start)) - start
    )),
    by = .(sample, TSR_FID)
  ]

  ## Calculate IQR.
  tss_position <- tss_position[order(sample, TSR_FID, tss_pos)]

  tss_position[,
    cum_sum := cumsum(score),
    by = .(sample, TSR_FID)
  ][,
    ecdf := ecdf_function(cum_sum),
    by = .(sample, TSR_FID)
  ]

  iqr_results <- tss_position[
    dplyr::between(ecdf, 0.1, 0.9)
  ][,
    .(iqr_min = min(cum_sum),
    iqr_max = max(cum_sum),
    iqr_width = max(tss_pos) - min(tss_pos),
    iqr_coords = str_c(seqnames, min(start), max(start), strand, sep = ":")),
    by = .(sample, TSR_FID)
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
#'
#' @rdname ecdf_function-function
#' @export

ecdf_function <- function(tss_vector) {
  ecdf_func <- ecdf(tss_vector)
  ecdf_values <- ecdf_func(tss_vector)
  return(ecdf_values)
}
