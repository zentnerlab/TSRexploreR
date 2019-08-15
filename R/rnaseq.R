
#' Add RNAseq Data
#'
#' Helper function to add RNA-seq data.
#'
#' @param experiment tsrexplorer object
#' @param rnaseq_count_matrix Raw counts in count matrix form from RNA-seq data
#'
#' @export
#' @rdname add_rna_seq-function

add_rna_seq <- function(experiment, rnaseq_count_matrix) {
	
	## Formatting count matrix for integration into tsrexplorer object.
	rnaseq_count_matrix <- as_tibble(rnaseq_count_matrix, .name_repair = "unique", rownames = "gene_id")

	## Adding RNA-seq count matrix to tsrexplorer object.
	experiment@raw_counts$RNAseq <- rnaseq_count_matrix

	return(experiment)
}

#' Add TSS Total Counts
#'
#' Add TSS feature counts data.
#'
#' @param experiment tsrexplorer object
#' @param tss_total_count_matrix Raw counts in count matrix form from TSS mapping data
#'
#' @export
#' @rdname add_tss_total-function

add_tss_total <- function(experiment, tss_total_count_matrix) {

	## Formatting count matrix for integration into tsrexplorer object.
	tss_total_count_matrix <- as_tibble(rnaseq_count_matrix, .name_repair = "unique", rownames = "gene_id")

	## Adding TSS total count matrix to tsrexplorer object.
	experiment@raw_counts$TSSs_total <- tss_total_count_matrix
}
