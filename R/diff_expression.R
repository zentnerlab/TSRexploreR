
#' edgeR Model for DE
#'
#' Find differential TSSs, TSRs, or features
#'
#' @import tibble
#' @importFrom edgeR DGEList filterByExpr calcNormFactors cpm estimateDisp glmQLFit
#' @importFrom dplyr select_at rename
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object after TMM normalization.
#' @param data_type Whether TSSs, TSRs, or feature counts should be analyzed.
#' @param samples Vector of sample names to analyze.
#' @param groups Vector of groups in correct factor notation.
#'
#' @return DGEList object with fitted model
#'
#' @rdname fit_edger_model-function
#'
#' @export

fit_edger_model <- function(experiment, data_type = c("tss", "tsr", "features"), samples = c(), groups = c()) {

	## Grab data from appropriate slot.
	sample_data <- experiment %>%
		extract_counts(data_type, samples) %>%
		bind_rows(.id = "sample") %>%
		as.data.table

	sample_data <- 

	## Select samples and turn to count matrix.
	selected_samples <- sample_data %>%
		select_at(c("position", samples)) %>%
		column_to_rownames("position") %>%
		as.matrix

	## Setting sample design.
	groups_factor <- factor(groups, levels = sort(unique(groups)))
	sample_design <- model.matrix(~ 0 + groups_factor)

	## Create DE fitted object.
	fitted_model <- selected_samples %>%
		DGEList(group = groups_factor) %>%
		.[filterByExpr(.), , keep.lib.sizes = FALSE] %>%
		calcNormFactors %>%
		estimateDisp(design = sample_design) %>%
		glmQLFit(design = sample_design)

	return(fitted_model)
}

#' Find Differential Expression
#'
#' Find differential TSSs, TSRs, or features from edgeR model.
#'
#' @import tibble
#' @importFrom edgeR glmQLFTest
#' @importFrom dplyr pull mutate
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#'
#' @param fit_edger_model edgeR differential expression model from fit_edger_model
#' @param data_type Whether the input was made from TSSs, TSRs, or features
#' @param compare_groups Vector of length two of the two groups to find differential TSRs
#'
#' @return tibble of differential TSRs
#'
#' @rdname differential_expression-function
#'
#' @export

differential_expression <- function(fit_edger_model, data_type = c("tss", "tsr", "feature"), compare_groups = c()) {
	
	## Set up contrasts.
	comparison_contrast <- fit_edger_model$samples %>%
		pull(group) %>%
		levels %>%
		as.numeric %>%
		length %>%
		numeric(length = .)

	comparison_contrast[compare_groups[1]] <- -1
	comparison_contrast[compare_groups[2]] <- 1

	## Differential expression
	diff_expression <- glmQLFTest(fit_edger_model, contrast = comparison_contrast)

	## Prepare tibble for export
	diff_expression <- diff_expression$table %>%
		as_tibble(.name_repair = "universal", rownames = "position") %>%
		mutate(FDR = p.adjust(PValue, method = "fdr")) %>%
		rename(log2FC = logFC)

	if (data_type %in% c("tss", "tsr")) {
		diff_expression <- separate(
			diff_expression,
			position,
			into = c("chr", "start", "end", "strand"),
			sep = "_"
		)
	} else if (data_type == "feature") {
		diff_expression <- rename(diff_expression, "gene_id" = position)
	}

	return(diff_expression)
}

#' Annotate Differential TSSs or TSRs
#'
#' Annotate Differential TSSs or TSRs to nearest gene or transcript.
#'
#' @import tibble
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom magrittr %>%
#'
#' @param differential_exp Tibble of differential TSSs or TSRs from differential_expression
#' @param annotation_file GTF genomic annotation file
#' @param feature_type Whether to annotate TSSs or TSRs relative to genes or transcripts
#' @param upstream Bases upstream of TSS
#' @param downstream Bases downstream of TSS
#'
#' @return Tibble of annotated differential TSSs or TSRs
#'
#' @rdname annotate_differential-function
#'
#' @export

annotate_differential <- function(
	differential_exp,
	annotation_file,
	feature_type = c("gene", "transcript"),
	upstream = 1000,
	downstream = 100
) {
	## Load genome annotation file as TxDb.
	genome_annotation <- makeTxDbFromGFF(annotation_file, "gtf")

	## Annotate differential TSRs.
	annotated_diff <- differential_exp %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
		annotatePeak(
			tssRegion = c(-upstream, downstream),
			TxDb = genome_annotation,
			sameStrand = TRUE,
			level = feature_type
		) %>%
		as_tibble(.name_repair = "universal")

	return(annotated_diff)
}

#' DE Volcano Plot
#'
#' Generate volcano plot for differential TSRs or Genes (RNA-seq)
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr case_when mutate
#'
#' @param differential_expression Tibble of differential TSRs or genes (RNA-seq) from differential_expression
#' @param log2fc_cutoff Log2 fold change cutoff for significance
#' @param fdr_cutoff FDR value cutoff for significance
#' @param ... Arguments passed to geom_point
#'
#' @return ggplot2 object of differential TSRs volcano plot.
#'
#' @rdname plot_volcano-function
#'
#' @export

plot_volcano <- function(
	differential_expression,
	log2fc_cutoff = 1,
	fdr_cutoff = 0.05,
	...
){

	## Annotate TSRs based on significance cutoff.
	diff_expression <- differential_expression %>%
		mutate(Change = case_when(
			log2FC >= log2fc_cutoff & FDR <= fdr_cutoff ~ "Increased",
			log2FC <= -log2fc_cutoff & FDR <= fdr_cutoff ~ "Decreased",
			TRUE ~ "Unchanged"
		)) %>%
		mutate(Change = factor(Change, levels = c("Decreased", "Unchanged", "Increased")))

	## Volcano plot of differential expression
	p <- ggplot(diff_expression, aes(x = log2FC, y = -log10(FDR))) +
		geom_point(aes(color = Change), ...) +
		scale_color_viridis_d() +
		theme_bw() +
		geom_vline(xintercept = -log2fc_cutoff, lty = 2) +
		geom_vline(xintercept = log2fc_cutoff, lty = 2) +
		geom_hline(yintercept = -log10(fdr_cutoff), lty = 2)

	return(p)
}

#' Export to clusterProfiler
#'
#' Export DEGs for use in clusterProfiler term enrichment.
#'
#' @import tibble
#' @importFrom dplyr select mutate case_when filter
#' 
#' @param annotated_de Annotated differential TSRs
#' @param log2fc_cutoff Log2 fold change cutoff for significance
#' @param fdr_cutoff FDR cutoff for significance
#'
#' @rdname export_for_enrichment-function
#'
#' @export

export_for_enrichment <- function(annotated_de, log2fc_cutoff = 1, fdr_cutoff = 0.05) {
	
	## Prepare data for export.
	export_data <- annotated_de %>%
		select(geneId, log2FC, FDR) %>%
		mutate(change = case_when(
			log2FC >= log2fc_cutoff & FDR <= fdr_cutoff ~ "increase",
			log2FC <= -log2fc_cutoff & FDR <= fdr_cutoff ~ "decrease",
			TRUE ~ "unchanged"
		)) %>%
		filter(change != "unchanged")

	return(export_data)
}
