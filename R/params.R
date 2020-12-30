
#' Function to hold common parameters.
#'
#' @param experiment TSRexploreR object.
#' @param genome_assembly Genome assembly in either a FASTA file or BSgenome object format.
#' @param genome_annotation Genome annotation in either a GTF, GFF3, or TxDb object format.
#' @param sample_sheet A sample sheet data.frame or tab delimited file.
#'   Must have the columns 'sample_name', 'file_1', and 'file_2'.
#'   Additional meta-data columns can be added with sample information such as condition and batch.
#' @param threshold TSSs or TSRs with a score below this value will be discarded.
#' @param samples A vector of sample names to analyze.
#' @param use_normalized Whether to use the normalized (TRUE) or raw (FALSE) counts.
#' @param ncol Integer specifying the number of columns to arrange multiple plots.
#' @param log2fc_cutoff Values below the |Log2 FC| cutoff will be discarded.
#' @param fdr_cutoff Values below the FDR cutoff will be discarded.
#' @param dominant If TRUE will discard any TSSs and/or TSRs that do not have the
#'   highest score per gene or transcript.
#' @param data_conditions Apply advanced conditions to the data.

common_params <- function(x) NULL
