source("setup.R")

test_that("Proper handling of TSS ranges.", {
  ## GRanges are stored to object.
  tsre <- tsr_explorer(TSSs["S288C_WT_1"])
  tsre@experiment$TSSs %>%
    expect_type("list") %>%
    expect_length(1) %>%
    expect_named("S288C_WT_1") %>%
    expect_identical(TSSs["S288C_WT_1"])

  ## GRanges are formatted to data.table.
  tsre <- format_counts(tsre, data_type="tss")
  tsre@counts$TSSs$raw %>%
    expect_length(1) %>%
    expect_named("S288C_WT_1")
  tsre@counts$TSSs$raw$S288C_WT_1 %>%
    expect_s3_class("data.table") %>%
    expect_named(c(
      "seqnames", "start", "end", "width",
      "strand", "score", "FHASH"
    ))
})

test_that("Proper handling of TSR ranges.", {
  tsre <- tsr_explorer(TSSs["S288C_WT_1"]) %>%
    format_counts(data_type="tss") %>%
    tss_clustering(threshold=3, max_distance=25)

  ## TSR GRanges are properly stored.
  tsre@experiment$TSRs %>%
    expect_type("list") %>%
    expect_length(1) %>%
    expect_named("S288C_WT_1")
  tsre@experiment$TSRs$S288C_WT_1 %>%
    expect_s4_class("GRanges")

  ## TSRs are converted to data.table.
  tsre@counts$TSRs$raw$S288C_WT_1 %>%
    expect_s3_class("data.table") %>%
    expect_named(c(
      "seqnames", "start", "end", "width",
      "strand", "score", "n_unique", "FHASH"
    ))
})

test_that("Association of TSRs with TSSs", {
  tsre <- tsr_explorer(TSSs["S288C_WT_1"]) %>%
    format_counts(data_type="tss") %>%
    tss_clustering(threshold=3, max_distance=25) %>%
    associate_with_tsr

  tsre@counts$TSSs$raw %>%
    expect_type("list") %>%
    expect_length(1) %>%
    expect_named("S288C_WT_1")
  tsre@counts$TSSs$raw$S288C_WT_1 %>%
    expect_s3_class("data.table") %>%
    expect_named(c(
      "seqnames", "start", "end", "width", "strand", "tsr_width", "tsr_score",
      "tsr_n_unique", "TSR_FHASH", "tsr_sample", "score", "FHASH"
    ))
})
