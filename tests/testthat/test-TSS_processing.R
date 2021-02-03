library("edgeR")
library("DESeq2")

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

test_that("Dominant TSS", {
  tsre <- TSSs["S288C_WT_1"] %>%
    tsr_explorer(genome_annotation=annotation) %>%
    format_counts(data_type="tss") %>%
    annotate_features(data_type="tss") %>%
    tss_clustering(threshold=3, max_distance=25) %>%
    associate_with_tsr %>%
    annotate_features(data_type="tsr")

  ## Dominant TSS per TSR.
  tss_tsr <- mark_dominant(tsre, data_type="tss")
  tss_tsr@counts$TSSs$raw %>%
    expect_length(1) %>%
    expect_type("list") %>%
    expect_named("S288C_WT_1")
  tss_tsr@counts$TSSs$raw$S288C_WT_1 %>%
    expect_s3_class("data.table") %>%
    {expect_true("dominant" %in% colnames(.))}

  ## Dominant TSS per gene.
  tss_gene <- mark_dominant(tsre, mark_per="gene", data_type="tss")
  tss_gene@counts$TSSs$raw %>%
    expect_length(1) %>%
    expect_type("list") %>%
    expect_named("S288C_WT_1")
  tss_gene@counts$TSSs$raw$S288C_WT_1 %>%
    expect_s3_class("data.table") %>%
    {expect_true("dominant" %in% colnames(.))}

  ## Dominant TSR per gene.
  tsr_gene <- mark_dominant(tsre, data_type="tsr")
  tsr_gene@counts$TSRs$raw %>%
    expect_length(1) %>%
    expect_type("list") %>%
    expect_named("S288C_WT_1")
  tsr_gene@counts$TSRs$raw$S288C_WT_1 %>%
    expect_s3_class("data.table") %>%
    {expect_true("dominant" %in% colnames(.))}

})
