#' Composition of Softclipped Bases
#'
#' @inheritParams common_params
#' @param n_bases Number of bases from -1 position to keep
#'
#' @export

softclip_composition <- function(
  experiment,
  samples="all",
  n_bases=NULL,
  ncol=3
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(is.null(n_bases) || is.count(n_bases))
  assert_that(is.count(ncol))

  ## Retrieve selected samples.
  if (all(samples == "all")) {
    select_samples <- experiment@experiment$TSSs
  } else {
    select_samples <- experiment@experiment$TSSs[samples]
  }

  ## Prepare data for softclip calculation.
  select_samples <- select_samples %>%
    map(as.data.table) %>%
    rbindlist(idcol="sample")

  select_samples[,
    c("seqnames", "start", "end", "strand", "width") := NULL
  ]

  ## Trim bases if n_bases is set.
  if (!is.null(n_bases)) {
    select_samples[
      n_soft > 0,
      seq_soft := str_sub(seq_soft, start=-n_bases, end=-1)
    ]
  }

  ## Get base content per soft-clip position.
  max_bases <- max(select_samples[, n_soft])

  select_samples[
    n_soft == 0,
    seq_soft := str_c(rep("-", max_bases), collapse="")
  ][
    n_soft > 0,
    seq_soft := str_pad(seq_soft, max_bases, "left", "-")
  ][,
    c(as.character(seq(-max_bases, -1, 1))) := tstrsplit(seq_soft, "")
  ][,
    c("rowid", "seq_soft", "n_soft") := list(.I, NULL, NULL)
  ]

  select_samples <- melt(
    select_samples, id.vars=c("sample", "rowid"),
    variable.name="position", value.name="base"
  )
  select_samples[,
    position := fct_relevel(position, as.character(seq(-1, -max_bases, -1)))
  ]

  ## Get table of base frequencies.
  select_samples <- select_samples[,
    .(count=.N),
    by=.(sample, position, base)
  ]

  ## Set sample order if required.
  if (!all(samples == "all")) {
    select_samples[, sample := factor(sample, levels=samples)]
  }

  ## Plot frequencies.
  p <- ggplot(select_samples, aes(x=.data$sample, y=.data$count)) +
    geom_col(aes(fill=.data$base), position="fill") +
    facet_wrap(~position, ncol=ncol)

  return(p)
}

#' Number of soft-clipped bases.
#'
#' @inheritParams common_params
#' @param n_bases Number of bases to plot
#'
#' @export

softclip_histogram <- function(
  experiment,
  samples="all",
  n_bases=NULL,
  ncol=3
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(is.null(n_bases) || is.count(n_bases))
  assert_that(is.count(ncol))

  ## Get samples.
  if (all(samples == "all")) {
    select_samples <- experiment@experiment$TSSs
  } else {
    select_samples <- experiment@experiment$TSSs[samples]
  }

  ## Prepare data for analysis.
  select_samples <- select_samples %>%
    map(as.data.table) %>%
    rbindlist(idcol="sample")
  select_samples[,
    c("seqnames", "start", "end", "strand", "width", "seq_soft") := NULL
  ]

  ## Set sample order if requested.
  if (!all(samples == "all")) {
    select_samples[, sample := factor(sample, levels=samples)]
  }

  ## Plot histogram.
  p <- ggplot(select_samples, aes(x=.data$n_soft)) +
    geom_histogram(aes(fill=.data$sample), binwidth=1) +
    facet_wrap(~ sample, ncol=ncol)

  return(p)
}
