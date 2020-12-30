context("Differential TSSs/TSRs")

source("setup.R")

sample_sheet <- data.frame(
  sample_name=c(sprintf("S288C_D_%s", seq_len(3)), sprintf("S288C_WT_%s", seq_len(3))),
  file_1=NA, file_2=NA,
  condition=c(rep("Diamide", 3), rep("Untreated", 3))
)

test_that("Differential TSSs", {
  tsre <- tsr_explorer(TSSs, sample_sheet=sample_sheet) %>%
    format_counts(data_type="tss")

  ## DEseq2.
  deseq2_columns <- c(
    "feature", "mean_expr", "log2FC", "stat", "pvalue", "padj",
    "seqnames", "start", "end", "strand"
  )
  deseq2_tsre <- fit_de_model(tsre, data_type="tss", formula= ~condition, method="DEseq2")
  expect_s4_class(deseq2_tsre@diff_features$TSSs$model, "DESeqDataSet")
  deseq2_tsre <- differential_expression(
    deseq2_tsre, data_type="tss",
    comparison_name="Diamide_vs_Untreated",
    comparison_type="name",
    comparison="condition_Untreated_vs_Diamide"
  )
  deseq2_tsre@diff_features$TSSs$results %>%
    expect_length(1) %>%
    expect_type("list") %>%
    expect_named("Diamide_vs_Untreated")
  deseq2_tsre@diff_features$TSSs$results$Diamide_vs_Untreated %>%
    expect_s3_class("data.table") %>%
    expect_named(deseq2_columns)

  ## edgeR.
  edger_columns <- c(
    "feature", "log2FC", "mean_expr", "pvalue", "padj", "seqnames", 
    "start", "end", "strand"
  )
  edger_tsre <- fit_de_model(tsre, data_type="tss", formula= ~condition, method="edgeR")
  expect_s4_class(edger_tsre@diff_features$TSSs$model, "DGEGLM")
  edger_tsre <- differential_expression(
    edger_tsre, data_type="tss",
    comparison_name="Diamide_vs_Untreated",
    comparison_type="contrast",
    comparison=c(1, -1)
  )
  edger_tsre@diff_features$TSSs$results$Diamide_vs_Untreated %>%
    expect_s3_class("data.table") %>%
    expect_named(edger_columns)
})
