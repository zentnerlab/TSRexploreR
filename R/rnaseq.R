
#' Add Feature Counts
#'
#' Helper function to add RNA-seq and 5' feature count data to tsrexplorer object.
#'
#' @param experiment tsrexplorer object
#' @param rnaseq_feature_counts Raw counts in count matrix form from RNA-seq data
#' @param five_prime_feature_counts Raw counts in count matrix form from 5' sequencing data
#'
#' @rdname add_feature_counts-function
#'
#' @export

add_feature_counts <- function(experiment, rnaseq_feature_counts, five_prime_feature_counts) {
  
  ## Formatting count matrix for integration into tsrexplorer object.
  rnaseq_feature_counts <- as_tibble(rnaseq_feature_counts, .name_repair = "unique", rownames = "gene_id")
  five_prime_feature_counts <- as_tibble(five_prime_feature_counts, .name_repair = "unique", rownames = "gene_id")

  ## Adding RNA-seq count matrix to tsrexplorer object.
  experiment@experiment$features <- list(
    "rna_seq" = rnaseq_feature_counts,
    "five_prime" = five_prime_feature_counts
  )

  return(experiment)
}
