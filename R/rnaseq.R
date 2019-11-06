
#' Add RNAseq Data
#'
#' Helper function to add RNA-seq data.
#'
#' @import tibble
#'
#' @param experiment tsrexplorer object
#' @param rnaseq_count_matrix Raw counts in count matrix form from RNA-seq data
#'
#' @rdname add_rnaseq_feature_counts-function
#'
#' @export

add_rnaseq_feature_counts <- function(experiment, rnaseq_feature_counts) {
	
	## Formatting count matrix for integration into tsrexplorer object.
	rnaseq_feature_counts <- as_tibble(rnaseq_feature_counts, .name_repair = "unique", rownames = "gene_id")

	## Adding RNA-seq count matrix to tsrexplorer object.
	experiment@raw_counts$RNAseq_features <- rnaseq_feature_counts

	return(experiment)
}

#' Add TSS Total Counts
#'
#' Add TSS feature counts data.
#'
#' @param experiment tsrexplorer object
#' @param tss_total_count_matrix Raw counts in count matrix form from TSS mapping data
#'
#' @rdname add_tss_feature_counts-function
#'
#' @export

add_tss_feature_counts <- function(experiment, tss_feature_counts) {

	## Formatting count matrix for integration into tsrexplorer object.
	tss_feature_counts <- as_tibble(tss_feature_counts, .name_repair = "unique", rownames = "gene_id")

	## Adding TSS total count matrix to tsrexplorer object.
	experiment@raw_counts$TSS_features <- tss_feature_counts

	return(experiment)
}
