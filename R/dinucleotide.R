#' Dinucleotide Analysis
#'
#' @description
#' Analysis of -1 and +1 dinucleotide frequencies.
#'
#' @include TSRexplore.R
#'
#' @inheritParams common_params
#' @param ... Arguments passed to geom_col
#'
#' @details
#' It has been shown in many organisms that particular base preferences exist at the
#'  -1 and +1 positions, where +1 is the TSS and -1 is the position immediately upstream.
#' This plotting function returns a ggplot2 barplot of -1 and +1 dinucleotide frequencies,
#'
#' 'genome_assembly' must be a valid genome assembly in either fasta or BSgenome format.
#' fasta formatted genome assemblies should have the file extension '.fasta' or '.fa'.
#' BSgenome assemblies are precompiled Bioconductor libraries for common organisms.
#'
#' A set of arguments to control data structure for plotting are included.
#' 'threshold' will define the minimum number of raw counts a TSS or TSR
#'  must have to be considered.
#' 'dominant' specifies whether only the dominant TSS should be considered 
#'   from the 'mark_dominant' function.
#' For TSSs this can be either dominant per TSR or gene.
#' 'data_conditions' allows for the advanced filtering, ordering, and grouping
#'   of data.
#'
#' @return ggplot2 object of dinucleotide plot.
#'   If 'return_table' is TRUE, a data.frame of underlying data is returned.
#'
#' @examples
#' data(TSSs)
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#'
#' tsre <- TSSs[1] %>%
#'   tsr_explorer(genome_annotation=annotation) %>%
#'   format_counts(data_type="tss")
#'
#' \donttest{plot_dinucleotide_frequencies(tsre)}
#'
#' @export

plot_dinucleotide_frequencies <- function(
  experiment,
  genome_assembly=NULL,
  samples="all",
  threshold=NULL,
  use_normalized=FALSE,
  dominant=FALSE,
  data_conditions=NULL,
  ncol=3,
  return_table=FALSE,
  ...
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is.null(genome_assembly) || is.character(genome_assembly) ||
    is(genome_assembly, "BSgenome")
  )
  assert_that(is.character(samples))
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.flag(use_normalized))
  assert_that(is.flag(dominant))
  assert_that(is.null(data_conditions) || is.list(data_conditions))
  assert_that(is.flag(return_table))

  ## Load genome assembly.
  fasta_assembly <- .prepare_assembly(genome_assembly, experiment)

  ## Get appropriate samples.
  select_samples <- experiment %>%
    extract_counts("tss", samples, use_normalized) %>%
    preliminary_filter(dominant, threshold)

  ## Apply conditions to data.
  select_samples <- condition_data(select_samples, data_conditions)

  ## Prepare samples for analysis.
  select_samples <- rbindlist(select_samples, idcol="sample")

  ## Extend ranges to capture base before TSS.
  select_samples <- select_samples %>%
    as_granges %>%
    anchor_3p %>%
    plyranges::mutate(width=2)

  ## Get dinucleotides.
  assembly_type <- case_when(
    is(fasta_assembly, "BSgenome") ~ "bsgenome",
    is(fasta_assembly, "FaFile") ~ "fafile"
  )

  seqs <- switch(
    assembly_type,
    "bsgenome"=BSgenome::getSeq(fasta_assembly, select_samples),
    "fafile"=Rsamtools::getSeq(fasta_assembly, select_samples)
  )

  seqs <- seqs %>%
    as.data.table %>%
    {cbind(as.data.table(select_samples), .)}
  setnames(seqs, old="x", new="dinucleotide")

  ## Find dinucleotide frequencies.
  data_grouping <- case_when(
    !is.null(data_conditions$quantiling) ~ "row_quantile",
    !is.null(data_conditions$grouping) ~ "row_groups",
    TRUE ~ "none"
  )

  freqs <- .calculate_freqs(seqs, data_grouping)

  ## Order dinucleotides by mean occurance.
  freqs[,
    mean_freq := mean(freqs), by=dinucleotide
  ][,
    rank := dense_rank(mean_freq)
  ][,
    dinucleotide := fct_reorder(dinucleotide, rank)
  ]

  ## Order samples if required.
  if (!all(samples == "all")) {
    freqs[, sample := factor(sample, levels=samples)]
  }

  ## Return table if requested.
  if (return_table) return(as.data.frame(freqs))

  ## PLot dinucleotide frequencies.
  p <- ggplot(freqs, aes(x=.data$dinucleotide, y=.data$freqs)) +
    geom_col(width=0.5, aes(fill=.data$freqs), ...) +
    theme_bw() +
    coord_flip() +
    labs(
      x="Dinucleotide",
      y="Frequency"
    )

  if (data_grouping != "none") {
    p <- p + facet_grid(as.formula(str_c(data_grouping, "~", "sample")))
  } else {
    p <- p + facet_wrap(~ sample, ncol=ncol)
  }

  return(p)

}

#' Internal DinucleotideFrequency Calculation
#'
#' @param seqs Dinucleotide sequences
#' @param grouping_status Whether there is quantiling or grouping.

.calculate_freqs <- function(seqs, grouping_status) {

  ## Get the dinucleotide frequencies.
  if (grouping_status == "none") {
    freqs <- seqs[,
      .(count=.N), by=.(sample, dinucleotide)
    ][,
      .(dinucleotide, count, freqs=count / sum(count)),
      by=sample
    ]
  } else {
    freqs <- seqs[,
      .(count=.N), by=c("sample", "dinucleotide", grouping_status)
    ][,
      .(dinucleotide, count, freqs=count / sum(count)),
      by=c("sample", grouping_status)
    ]
  }

  return(freqs)
}
