
#' Number of Softclipped Bases
#'
#' @param experiment tsr explorer object
#'
#' @export

softclip_histogram <- function(
  experiment,
  samples="all",
  n_bases=NULL
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))

  ## Retrieve selected samples.
  if (samples == "all") {
    select_samples <- experiment@experiment$TSSs
  } else {
    select_samples <- experiment@experiment$TSSs[samples]
  }

  ## Prepare data for softclip calculation.
  select_samples <- select_samples %>%
    map(as.data.table) %>%
    rbindlist(idcol="sample")

  ## Trim bases if n_bases is set.
  if (!is.null(n_bases)) {
    select_samples[
      n_soft > 0,
      seq_soft := str_sub(seq_soft, start=-n_bases, end=-1)
    ]
  }

  ## Get table of base frequencies.
  freqs <- select_samples[, .(count=.N), by=.(sample, seq_soft)]
  freqs[,
    freq := count/sum(count), by=sample
  ][,
    seq_soft := fct_reorder(seq_soft, count, .desc=TRUE)
  ]

  kable(freqs, format="simple", caption="Soft-clipped Base Frequency")

  ## Plot frequencies.
  p <- ggplot(freqs, aes(x=sample, y=count)) +
    geom_col(aes(fill=seq_soft), position="fill") +
    coord_flip()

  return(p)
}
