
#' Add RNAseq Data
#'
#' Helper function to add RNA-seq data.
#'
#' @import tibble
#'
#' @param experiment tsrexplorer object
#' @param rnaseq_count_matrix Raw counts in count matrix form from RNA-seq data
#'
#' @export
#' @rdname add_rnaseq_feature_counts-function

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
#' @export
#' @rdname add_tss_feature_counts-function

add_tss_feature_counts <- function(experiment, tss_feature_counts) {

	## Formatting count matrix for integration into tsrexplorer object.
	tss_feature_counts <- as_tibble(tss_feature_counts, .name_repair = "unique", rownames = "gene_id")

	## Adding TSS total count matrix to tsrexplorer object.
	experiment@raw_counts$TSS_features <- tss_feature_counts

	return(experiment)
}

#' RNA-seq TSSs total Correlation Matrix
#'
#' Correlation between total TSSs and RNA-seq
#'
#' @import tibble
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom dplyr left_join
#'
#' @param experiment tsrexplorer object with RNA-seq and TSS total added
#' @param corr_metric Correlation metrix to use ("pearon", "spearman")
#'
#' @return ggplot2 object of RNA-seq versus TSS mapping.
#'
#' @export
#' @rdname plot_rnaseq_corr-function

plot_rnaseq_corr <- function(experiment, corr_metrix = c("pearson", "spearman")) {
	
	## Preparing data for plotting.
	rnaseq_corr <- left_join(experiment@raw_counts$RNAseq, experiment@raw_counts$TSSs_total, by = "gene_id") %>%
		column_to_rownames("gene_id") %>%
		as.matrix %>%
		DGEList %>%
		calcNormFactors %>%
		cpm %>%
		cor(method = corr_metric) %>%
		as_tibble(rownames = "sample_1", .name_repair="unique") %>%
		gather(key = "sample_2", value = corr_metric, -sample_1) %>%
		mutate(corr_metric = round(corr_metric, 3))

	## Plotting correlation.
	p <- ggplot(rnaseq_corr, aes(x=sample_1, y=sample_2, fill=corr_metric, label=corr_metric)) +
		geom_tile(color="white", lwd=0.5) +
		geom_label(color="white", label.size=NA, fill=NA) +
		scale_fill_viridis_c(limits = c(0, 1), name=corr_metric, direction = -1) +
		theme_minimal() +
		theme(
			axis.text.x=element_text(angle=45, hjust=1),
			panel.grid=element_blank(),
			axis.title=element_blank()
		)

	return(p)
}
