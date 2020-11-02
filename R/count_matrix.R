
#' Generate Count Matrix.
#'
#' @param count_data Names list of counts.
#' @param data_type Currently supports 'tss' only.
#' @param use_normalized Whether to use normalized counts.

.count_matrix <- function(
  count_data,
  data_type=c("tss", "tsr"),
  use_normalized=FALSE
) {

  ## Input check.
  assert_that(is.list(count_data) && has_attr(count_data, "names"))
  data_type <- match.arg(str_to_lower(data_type), "tss")
  assert_that(is.flag(use_normalized))

  ## Change scores to normalized scores if requested.
  if (use_normalized) {
    walk(count_data, ~.x[, score := normalized_score])
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
