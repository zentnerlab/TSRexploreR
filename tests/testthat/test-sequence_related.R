context("Sequence Retrieval Functions")
source("setup.R")

test_that("TSS sequence colormap and sequence logo", {
  ## Sequence DataFrame.
  tsre <- tsr_explorer(TSSs["S288C_WT_1"], genome_assembly=assembly) %>%
    format_counts(data_type="tss")

  ## Sequence colormap.
  tsre %>%
    plot_sequence_colormap %>%
    expect_s3_class("ggplot")

  ## Sequence logo.
  tsre %>%
    plot_sequence_logo %>%
    expect_s3_class("ggplot")
})

test_that("Dinucleotide frequencies", {
  ## DataFrame with frequencies.
  tsre <- tsr_explorer(TSSs["S288C_WT_1"], genome_assembly=assembly) %>%
    format_counts(data_type="tss")

  ## Dinucleotide plot.
  tsre %>%
    plot_dinucleotide_frequencies %>%
    expect_s3_class("ggplot")
})
