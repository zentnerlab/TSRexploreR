
#' Generate Count Matrix.
#'
#' @inheritParams common_params
#' @param count_data Named list of counts.
#' @param data_type Whether to get matrix of TSS ('tss') or TSR ('tsr') counts.

.count_matrix <- function(
  count_data,
  data_type=c("tss", "tsr"),
  use_normalized=FALSE,
  threshold=NULL
) {

  ## Input check.
  assert_that(is.list(count_data) && has_attr(count_data, "names"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr")
  )
  assert_that(is.flag(use_normalized))
  assert_that(
    is.null(threshold) ||
    (is.numeric(threshold) && threshold > 0)
  )

  ## Change scores to normalized scores if requested.
  if (use_normalized) {
    walk(count_data, ~.x[, score := normalized_score])
  }

  ## If making a TSR matrix merge overlapping ranges.
  if (data_type == "tsr") {

    # Reduce overlapping ranges.
    merged_ranges <- count_data %>%
      map(as_granges) %>%
      bind_ranges %>%
      reduce_ranges_directed %>%
      as.data.table

    # Prepare reduced ranges for overlap with original ranges.
    setkey(merged_ranges, seqnames, strand, start, end)
    merged_ranges[,
      FHASH := str_c(seqnames, start, end, strand, sep=":"),
      by=seq_len(nrow(merged_ranges))
    ][,
      width := NULL
    ]

    # Get the aggregated score for reduced ranges.
    count_data <- map(count_data, function(x) {
      setkey(x, seqnames, strand, start, end)
      x <- foverlaps(x, merged_ranges)
      x <- x[, .(score=sum(score)), by=FHASH]
      return(x)
    })
    
  }

  ## Create the count matrix.
  count_data <- rbindlist(count_data, idcol="sample")[,
    .(FHASH, sample, score)
  ]
  count_data <- dcast(count_data, FHASH ~ sample, value.var="score", fill=0)

  count_data <- count_data %>%
    column_to_rownames("FHASH") %>%
    as.matrix

  ## Return the count matrix.
  return(count_data)
}
