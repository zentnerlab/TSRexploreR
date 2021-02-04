#' Composition of Softclipped Bases
#'
#' @description
#' Stacked barplot of base composition of softclipped bases
#'   upstream and adjacent to the TSS.
#'
#' @inheritParams common_params
#' @param n_bases Number of bases from -1 position to keep
#' @param ... Arguments passed to geom_col
#'
#' @details
#' Cap trapping and TSRT based 5' mapping methods have shown
#'   a preponderance of softclipped Gs immediatly upstream and
#'   adjacent to TSSs.
#' The source of this extra G is hypothesized to be reverse
#'   transcription of the 5' cap to a C.
#' In addition to this extra base, TSRT based methods have been
#'   shown to add in upwards of 2-4 additional bases.
#'
#' This function generates a stacked barplot for each position
#'   upstream of the TSS that gives the relative base composition of
#'   each softclipped base if present.
#' 'n_bases' determines how far upstream of the TSS is checked for
#'   sofclipped bases.
#'
#' If 'return_table' is TRUE, a data.frame is returned with the
#'   underlying numbers used in the plot.
#'
#' @return A ggplot2 stacked barplot of softclip base composition.
#'   If 'return_table' is TRUE a data.frame is returned instead.
#'
#' @seealso
#' \code{\link{import_bams}} to import BAMs.
#' \code{\link{G_correction}} to correct for incidentally templated 
#'   spurrious 5' Gs.
#'
#' @examples
#' bam_file <- system.file("extdata", "S288C.bam", package="TSRexploreR")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#' samples <- data.frame(sample_name="S288C", file_1=bam_file, file_2=NA)
#'
#' tsre <- tsr_explorer(sample_sheet=samples, genome_assembly=assembly) %>%
#'   import_bams(paired=TRUE)
#' p <- softclip_composition(tsre)
#'
#' @export

softclip_composition <- function(
  experiment,
  samples="all",
  n_bases=NULL,
  ncol=3,
  return_table=FALSE,
  ...
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(is.null(n_bases) || is.count(n_bases))
  assert_that(is.count(ncol))
  assert_that(is.flag(return_table))

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

  ## Return table if requested.
  if (return_table) return(as.data.frame(select_samples))

  ## Plot frequencies.
  p <- ggplot(select_samples, aes(x=.data$sample, y=.data$count)) +
    geom_col(aes(fill=.data$base), position="fill") +
    facet_wrap(~position, ncol=ncol)

  return(p)
}

#' Number of soft-clipped bases.
#'
#' @description
#' Histogram of the number of softlcipped bases upstream and adjacent to
#'   the TSS.
#'
#' @inheritParams common_params
#' @param n_bases Number of bases to plot
#'
#' @details
#' Cap trapping and TSRT based 5' mapping methods have shown
#'   a preponderance of softclipped Gs immediatly upstream and
#'   adjacent to TSSs.
#' The source of this extra G is hypothesized to be reverse
#'   transcription of the 5' cap to a C.
#' In addition to this extra base, TSRT based methods have been
#'   shown to add in upwards of 2-4 additional bases.
#'
#' This function creates a histogram for each position upstream of
#'   the TSS, up to 'n_bases' away.
#'
#' @return ggplot2 histogram of softclipped base numbers.
#'
#' @seealso
#' \code{\link{import_bams}} to import BAMs.
#' \code{\link{G_correction}} to correct for incidentally templated 
#'   spurrious 5' Gs.
#'
#' @examples
#' bam_file <- system.file("extdata", "S288C.bam", package="TSRexploreR")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#' samples <- data.frame(sample_name="S288C", file_1=bam_file, file_2=NA)
#'
#' tsre <- tsr_explorer(sample_sheet=samples, genome_assembly=assembly) %>%
#'   import_bams(paired=TRUE)
#' p <- softclip_histogram(tsre)
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
