context("Genome Annotation and Assembly")

source("setup.R")
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library("BSgenome.Scerevisiae.UCSC.sacCer3")

test_that("Import of assemblies", {
  ## FASTA assembly.
  fasta_tsre <- tsr_explorer(genome_assembly=assembly)
  expect_s4_class(fasta_tsre@meta_data$genome_assembly, "FaFile")

  ## BSgenome assembly.
  bsgenome_tsre <- tsr_explorer(genome_assembly=BSgenome.Scerevisiae.UCSC.sacCer3)
  expect_s4_class(bsgenome_tsre@meta_data$genome_assembly, "BSgenome")
})

test_that("Import of annotations", {
  ## GTF annotation.
  gtf_tsre <- tsr_explorer(genome_annotation=annotation)
  expect_s4_class(gtf_tsre@meta_data$genome_annotation, "TxDb")

  ## TxDb annotation.
  txdb_tsre <- tsr_explorer(genome_annotation=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
  expect_s4_class(txdb_tsre@meta_data$genome_annotation, "TxDb")
})

test_that("TSS annotation", {
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
