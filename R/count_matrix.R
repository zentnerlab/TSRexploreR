
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
    walk(select_samples, ~.x[, score := normalized_score])
  }

  ## Create the count matrix.
  select_samples <- rbindlist(select_samples, idcol="sample")[
    .(FHASH, sample, score)
  ]
  select_samples <- dcast(select_samples, FHASH ~ sample, value.var="score")
  setnafill(select_samples, fill=0)

  select_samples <- select_samples %>%
    column_to_rownames("FHASH") %>%
    as.matrix

  ## Return the count matrix.
  return(select_samples)
}
