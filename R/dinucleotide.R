
#' Dinucleotide Analysis
#'
#' @description
#' Analysis of -1 and +1 dinucleotide frequencies.
#'
#' @include tsrexplorer.R
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param samples Either 'all' to plot all samples or a vector of sample names
#' @param genome_assembly fasta file or BSgenome of genome assembly
#' @param threshold Raw TSS count threshold
#' @param dominant Consider only dominant TSSs
#' @param data_conditions Condition the data
#'
#' @details
#' It has been shown in many organisms that particular base preferences exist at the
#'  -1 and +1 positions, where +1 is the TSS and -1 is the position immediately upstream.
#' This function returns the dinucleotide frequencies at each TSS.
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
#' For TSSs this can be either dominant per TSR or gene, and for TSRs
#'   it is just the dominant TSR per gene. (qq should this just be TSSs?)
#' 'data_conditions' allows for the advanced filtering, ordering, and grouping
#'   of data.
#'
#' @return DataFrame with dinucleotide frequencies
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "tsrexplorer")
#' freqs <- dinucleotide_frequencies(tsre_exp, genome_assembly = assembly)
#'
#' @seealso
#' \code{\link{plot_dinculeotide_frequencies}} to plot the dinucleotide frequencies.
#'
#' @rdname dinucleotide_frequencies-function
#' @export

dinucleotide_frequencies <- function(
  experiment,
  genome_assembly,
  samples = "all",
  threshold = 1,
  dominant = FALSE,
  data_conditions = NA
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(genome_assembly) | is(genome_assembly, "BSgenome"))
  assert_that(is.character(samples))
  assert_that(is.count(threshold) && threshold > 0)
  assert_that(is.flag(dominant))
  if (all(!is.na(data_conditions)) && !is(data_conditions, "list")) {
    stop("data_conditions must be a list of values")
  }

  ## Load genome assembly.
  assembly_type <- case_when(
    is(genome_assembly, "character") ~ "character",
    is(genome_assembly, "BSgenome") ~ "bsgenome"
  )

  fasta_assembly <- switch(
    assembly_type,
    "character"=FaFile(genome_assembly),
    "bsgenome"=genome_assembly
  )

  ## Get appropriate samples.
  select_samples <- extract_counts(experiment, "tss", samples)

  ## Preliminary filtering of data.
  select_samples <- preliminary_filter(select_samples, dominant, threshold)

  ## Apply conditions to data.
  if (all(!is.na(data_conditions))) {
    select_samples <- do.call(group_data, c(list(signal_data = select_samples), data_conditions))
  }

  ## Prepare samples for analysis.
  select_samples <- rbindlist(select_samples, idcol = "sample")
  select_samples <- select_samples[, .(seqnames, start, end, strand, sample)]

  ## Extend ranges to capture base before TSS.
  select_samples <- select_samples %>%
    as_granges %>%
    anchor_3p %>%
    plyranges::mutate(width=2)

  ## Get dinucleotides.
  seqs <- select_samples %>%
    {getSeq(fasta_assembly, .)} %>%
    as.data.table %>%
    {cbind(as.data.table(select_samples), .)}
  setnames(seqs, old = "x", new = "dinucleotide")

  ## Find dinucleotide frequencies.
  groupings <- any(names(data_conditions) %in% c("quantile_by", "grouping"))

  if (!groupings) {
    freqs <- seqs[,
      .(count = .N), by = .(sample, dinucleotide)
    ][,
      .(dinucleotide, count, freqs = count / sum(count)),
      by = sample
    ]
  } else {
    freqs <- seqs[,
      .(count = .N), by = .(sample, dinucleotide, grouping)
    ][,
      .(dinucleotide, count, freqs = count / sum(count)),
      by = .(sample, grouping)
    ]
  }

  ## Order samples if required.
  if (!all(samples == "all")) {
    freqs[, sample := factor(sample, levels = samples)]
  }

  ## Prepare DataFrame to return.
  freqs_df <- DataFrame(freqs)
  metadata(freqs_df)$groupings <- groupings
  metadata(freqs_df)$threshold <- threshold
  
  return(freqs_df)
}

#' Plot Dinucleotide Frequencies
#'
#' @description
#' Plot results from dinucleotide analysis
#'
#' @param dinucleotide_frequencies tibble from dinucleotide_frequencies analysis
#' @param ncol Number of columns to plot if quantiles are not specified
#' @param ... Arguments passed to geom_col
#'
#' @details
#' This plotting function returns a ggplot2 barplot of -1 and +1 dinucleotide frequencies, 
#'   where +1 is the TSS and -1 is the position immediately upstream.
#' The results of the 'dinucleotide_frequencies' function are used as input.
#'
#' @return ggplot2 object of dinucleotide frequencies plot
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "tsrexplorer")
#' freqs <- dinucleotide_frequencies(tsre_exp, genome_assembly = assembly)
#' plot_dinucleotide_frequencies(freqs)
#'
#' @seealso
#' \code{\link{dinucleotide_frequencies}} to calculate the dinucleotide frequencies.
#'
#' @rdname plot_dinucleotide_frequencies-function
#' @export

plot_dinucleotide_frequencies <- function(
  dinucleotide_frequencies,
  ncol = 1,
  ...
) {

  ## Check inputs.
  assert_that(is(dinucleotide_frequencies, "DataFrame"))
  assert_that(is.count(ncol) && ncol > 0)

  ## Pull out some info from the dataframe.
  groupings <- metadata(dinucleotide_frequencies)$groupings

  ## Convert dataframe to data.table.
  freqs <- as.data.table(dinucleotide_frequencies)
  
  ## Set factor order for dinucleotides.
  if (!groupings) {
    freqs <- freqs[,
      .(sample, count, freqs, mean_freqs = mean(freqs)),
      by = dinucleotide
    ]
  } else {
    freqs <- freqs[,
      .(sample, count, freqs, grouping, mean_freqs = mean(freqs)),
      by = dinucleotide
    ]
  }

  freqs[,
    rank := dense_rank(mean_freqs)
  ][,
    dinucleotide := fct_reorder(factor(dinucleotide), rank)
  ][,
    c("mean_freqs", "rank") := NULL
  ]

  ## Plot dinucleotide frequencies.

  p <- ggplot(freqs, aes(x = .data$dinucleotide, y = .data$freqs)) +
    geom_col(width = 0.5, aes(fill = .data$freqs), ...) +
    theme_bw() +
    coord_flip() +
    labs(
      x="Dinucleotide",
      y="Frequency"
    )

  if (groupings) {
    p <- p + facet_wrap(grouping ~ sample)
  } else {
    p <- p + facet_wrap(~ sample, ncol = ncol)
  }

  return(p)
}
