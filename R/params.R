
#' Function to hold common parameters.
#'
#' @param experiment TSRexploreR object.
#' @param genome_assembly Genome assembly in FASTA or BSgenome format.
#' @param genome_annotation Genome annotation in GTF, GFF3, or TxDb format.
#' @param sample_sheet A sample sheet data.frame or tab delimited file.
#'   Must have the columns 'sample_name', 'file_1', and 'file_2'.
#'   Additional meta-data columns can be added with sample information such as condition and batch.
#' @param threshold TSSs or TSRs with a score below this value will not be considered.
#' @param samples A vector of sample names to analyze.
#' @param use_normalized Whether to use the normalized (TRUE) or raw (FALSE) counts.
#' @param ncol Integer specifying the number of columns to arrange multiple plots.
#' @param log2fc_cutoff Differential features not meeting this log2(fold change) threshold will not be considered.
#' @param fdr_cutoff Differential features not meeting this significance threshold will not be considered.
#' @param dominant If TRUE, will only consider the highest-scoring TSS per gene, transcript, or TSR or 
#'   highest-scoring TSR per gene or transcript.
#' @param exclude_antisense Remove antisense TSSs/TSRs prior to analysis.
#' @param data_conditions Apply advanced conditions to the data.
#' @param rasterize Rasterize a ggplot.
#' @param raster_dpi If rasterization is set, this controls the rasterization DPI.
#' @param return_table Return a table of results instead of a plot.

common_params <- function(
  experiment,
  genome_assembly,
  genome_annotation,
  sample_sheet,
  threshold,
  samples,
  use_normalized,
  ncol,
  log2fc_cutoff,
  fdr_cutoff,
  dominant,
  exclude_antisense,
  data_conditions,
  rasterize,
  raster_dpi,
  return_table
) {
NULL
}
