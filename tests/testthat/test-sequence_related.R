context("Sequence Retrieval Functions")
source("setup.R")

test_that("TSS sequence colormap and sequence logo", {
  ## Sequence DataFrame.
  seqs <- tsr_explorer(TSSs["S288C_WT_1"], genome_assembly=assembly) %>%
    format_counts(data_type="tss") %>%
    tss_sequences
  seqs %>%
    expect_s4_class("DataFrame") %>%
    expect_named(c("sample", "FHASH", "plot_order", "tss", "sequence", "score"))

  ## Sequence colormap.
  seqs %>%
    plot_sequence_colormap %>%
    expect_s3_class("ggplot")

  ## Sequence logo.
  seqs %>%
    plot_sequence_logo %>%
    expect_s3_class("ggplot")
})

test_that("Dinucleotide frequencies", {
  ## DataFrame with frequencies.
  seqs <- tsr_explorer(TSSs["S288C_WT_1"], genome_assembly=assembly) %>%
    format_counts(data_type="tss") %>%
    dinucleotide_frequencies
  seqs %>%
    expect_s4_class("DataFrame") %>%
    expect_named(c("sample", "dinucleotide", "count", "freqs"))

  ## Dinucleotide plot.
  seqs %>%
    plot_dinucleotide_frequencies %>%
    expect_s3_class("ggplot")
})
