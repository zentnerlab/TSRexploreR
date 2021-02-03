library("magrittr")

bam_file <- system.file("extdata", "S288C.bam", package="TSRexploreR")
samples <- data.frame(sample_name="S288C", file_1=bam_file, file_2=NA)

test_that("Import and Process BAM", {
  tsre <- tsr_explorer(sample_sheet=samples, genome_assembly=assembly)

  ## Bam import.
  tsre <- import_bams(tsre, paired=TRUE)
  tsre@experiment$TSSs %>%
    expect_length(1) %>%
    expect_type("list") %>%
    expect_named("S288C")
  tsre@experiment$TSSs$S288C %>%
    expect_s4_class("GRanges") %>%
    {expect_named(mcols(.), c("seq_soft", "n_soft"))}

  ## G correction.
  tsre <- G_correction(tsre)
  tsre@experiment$TSSs %>%
    expect_length(1) %>%
    expect_type("list") %>%
    expect_named("S288C")
  tsre@experiment$TSSs$S288C %>%
    expect_s4_class("GRanges") %>%
    {expect_named(mcols(.), c("seq_soft", "n_soft"))}

  ## TSS aggregation.
  tsre <- tss_aggregate(tsre)
  tsre@experiment$TSSs %>%
    expect_length(1) %>%
    expect_type("list") %>%
    expect_named("S288C")
  tsre@experiment$TSSs$S288C %>%
    expect_s4_class("GRanges") %>%
    {expect_named(mcols(.), "score")}
})
