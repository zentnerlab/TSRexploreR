context("Genome Annotation and Assembly")

library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library("BSgenome.Scerevisiae.UCSC.sacCer3")
assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")

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
