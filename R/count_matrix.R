
#' Generate Count Matrix.
#'
#' @param count_data Names list of counts.
#' @param use_normalized Whether to use normalized counts.

.count_matrix <- function(
  count_data,
  data_type=c("tss", "tsr"),
  use_normalized=FALSE
) {

  ## Input check.
  assert_that(is.list(count_data) && has_attr(count_data, "names"))
  assert_that(is.flag(use_normalized))

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
      FHASH := str_c(seqnames, start, end, strand),
      by=seq_len(nrow(merged_ranges))
    ][,
      width := NULL
    ]

    # Get the aggregated score for reduced ranges.
    count_data <- map(count_data, function(x) {
      setkey(x, seqnames, strand, start, end)
      x <- foverlaps(x, merged_ranges)
      x[,
        c("FHASH", "i.start", "i.end") := list(
          i.FHASH, NULL, NULL
        )
      ]
      x <- x[, .(score=sum(score)), by=FHASH]
      return(x)
    })
    
  }

  ## Create the count matrix.
  count_data <- rbindlist(count_data, idcol="sample")[,
    .(FHASH, sample, score)
  ]
  count_data <- dcast(count_data, FHASH ~ sample, value.var="score")

  fill_cols <- colnames(count_data)[!colnames(count_data) == "FHASH"]
  setnafill(count_data, fill=0, cols=fill_cols)

  count_data <- count_data %>%
    column_to_rownames("FHASH") %>%
    as.matrix

  ## Return the count matrix.
  return(count_data)
}
