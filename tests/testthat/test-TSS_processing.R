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

test_that("TSS annotation", {
  annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
  tsre <- tsr_explorer(TSSs["S288C_WT_1"], genome_annotation=annotation) %>%
    format_counts(data_type="tss")

  expected_columns <- c(
    "seqnames", "start", "end", "width", "strand", "score", "FHASH", 
    "annotation", "geneChr", "geneStart", "geneEnd", "geneLength", 
    "geneStrand", "geneId", "distanceToTSS", "simple_annotations"
  )

  ## Annotate based on gene.
  gene_columns <- c(
    "seqnames", "start", "end", "width", "strand", "score", "FHASH",
    "annotation", "geneChr", "geneStart", "geneEnd", "geneLength",
    "geneStrand", "geneId", "distanceToTSS", "simple_annotations"
  )
  gene_anno <- annotate_features(tsre, data_type="tss", feature_type="gene")
  gene_anno@counts$TSSs$raw %>%
    expect_type("list") %>%
    expect_length(1) %>%
    expect_named("S288C_WT_1")
  gene_anno@counts$TSSs$raw$S288C_WT_1 %>%
    expect_s3_class("data.table") %>%
    expect_named(gene_columns)

  ## Annotate based on transcript.
  tx_columns <- c(
    "seqnames", "start", "end", "width", "strand", "score", "FHASH", 
    "annotation", "geneChr", "geneStart", "geneEnd", "geneLength", 
    "geneStrand", "geneId", "transcriptId", "distanceToTSS", "simple_annotations"
  )
  tx_anno <- annotate_features(tsre, data_type="tss", feature_type="transcript")
  tx_anno@counts$TSSs$raw %>%
    expect_type("list") %>%
    expect_length(1) %>%
    expect_named("S288C_WT_1")
  tx_anno@counts$TSSs$raw$S288C_WT_1 %>%
    expect_s3_class("data.table") %>%
    expect_named(tx_columns)
})
