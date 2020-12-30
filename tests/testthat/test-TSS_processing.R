context("TSS Processing")
source("setup.R")

test_that("TSS Normalization", {
  sample_sheet <- data.frame(
    sample_name=c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3"),
    file_1=NA, file_2=NA
  )
  tsre <- tsr_explorer(
    TSSs[c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3")],
    sample_sheet=sample_sheet
  ) %>%
    format_counts(data_type="tss")

  ## CPM normalization.
  cpm_norm <- normalize_counts(tsre, data_type="tss", method="CPM")
  cpm_norm@counts$TSSs$raw %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3")) %>%
    walk(expect_s3_class, "data.table") %>%
    walk(expect_named, c(
      "seqnames", "start", "end", "width", "strand",
      "FHASH", "normalized_score", "score"
    ))

  ## edgeR normalization.
  edger_norm <- normalize_counts(tsre, data_type="tss", method="edgeR")
  edger_norm@counts$TSSs$raw %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3")) %>%
    walk(expect_s3_class, "data.table") %>%
    walk(expect_named, c(
      "seqnames", "start", "end", "width", "strand",
      "FHASH", "normalized_score", "score"
    ))

  ## DESeq2 normalization.
  deseq2_norm <- normalize_counts(tsre, data_type="tss", method="DESeq2")
  deseq2_norm@counts$TSSs$raw %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3")) %>%
    walk(expect_s3_class, "data.table") %>%
    walk(expect_named, c(
      "seqnames", "start", "end", "width", "strand",
      "FHASH", "normalized_score", "score"
    ))
})
