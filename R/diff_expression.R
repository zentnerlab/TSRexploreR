
#' edgeR Model for DE
#'
#' Find differential TSSs, TSRs, or features
#'
#' @import tibble
#' @importFrom edgeR DGEList filterByExpr calcNormFactors cpm estimateDisp glmQLFit
#' @importFrom dplyr select_at rename
#' @importFrom magrittr %>%
#' @importFrom forcats fct_inorder
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

fit_edger_model <- function(
	experiment, data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	samples, groups
) {
	## Design table.
	design <- data.table("samples" = samples, "groups" = groups)
	design[, groups := fct_inorder(as.character(groups))]

	## Grab data from appropriate slot.
	sample_data <- extract_matrix(experiment, data_type, design[["samples"]])

	## Filter out features with low counts.
	sample_data <- sample_data[filterByExpr(assay(sample_data, "counts")), ]

	## Setting sample design.
	sample_design <- model.matrix(~ 0 + design[["groups"]])

	## Create DE fitted object.
	fitted_model <- assay(sample_data, "counts") %>%
		set_rownames(rowRanges(sample_data)$FID) %>%
		DGEList(group = design[["groups"]]) %>%
		calcNormFactors %>%
		estimateDisp(design = sample_design) %>%
		glmQLFit(design = sample_design)

	## Store model in tsrexplorer object.
	if (data_type == "tss") {
		experiment@diff_features$TSSs$model <- fitted_model
		experiment@diff_features$TSSs$design <- design
	} else if (data_type == "tsr") {
		experiment@diff_features$TSRs$model <- fitted_model
		experiment@diff_features$TSRs$design <- design
	} else if (data_type == "tss_features") {
		experiment@diff_features$TSS_features$model <- fitted_model
		experiment@diff_features$TSS_features$design <- design
	} else if (data_type == "tsr_features") {
		experiment@diff_features$TSR_features$model <- fitted_model
		experiment@diff_features$TSR_features$design <- design
	}

	return(experiment)
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
#' @importFrom purrr map_dbl
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

differential_expression <- function(
	experiment, data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	compare_groups, fdr_cutoff = 0.05, log2fc_cutoff = 1) {
	
	## Grab appropriate model.
	if (data_type == "tss") {
		edger_model <- experiment@diff_features$TSSs$model
		edger_design <- experiment@diff_features$TSSs$design
		anno_data <- as.data.table(rowRanges(experiment@counts$TSSs$matrix))
	} else if (data_type == "tsr") {
		edger_model <- experiment@diff_features$TSRs$model
		edger_design <- experiment@diff_features$TSRs$design
		anno_data <- as.data.table(rowRanges(experiment@counts$TSRs$matrix))
	} else if (data_type == "tss_features") {
		edger_model <- experiment@diff_features$TSS_features$model
		edger_design <- experiment@diff_features$TSS_features$design
		anno_data <- as.data.table(rowData(experiment@counts$TSS_features$matrix))
	} else if (data_type == "tsr_features") {
		edger_model <- experiment@diff_features$TSR_features$model
		edger_design <- experiment@diff_features$TSR_features$design
		anno_data <- as.data.table(rowData(experiment@counts$TSR_features$matrix))
	}

	## Set up contrasts.
	comparison_levels <- edger_design %>%
		pull(groups) %>%
		levels

	comparison_contrast <- comparison_levels %>%
		map_dbl(function(x) {
			if (x == compare_groups[1]) {
				return_value <- 1
			} else if (x == compare_groups[2]) {
				return_value <- -1
			} else {
				return_value <- 0
			}
			return(return_value)
		})

	## Differential expression
	diff_expression <- glmQLFTest(edger_model, contrast = comparison_contrast)

	diff_expression <- diff_expression$table %>%
		rownames_to_column("FID") %>%
		as.data.table

	setnames(diff_expression, old = "logFC", new = "log2FC")

	comparison_name <- str_c(compare_groups[1], "_vs_", compare_groups[2])
	diff_expression[, c("FDR", "sample") := list(
		p.adjust(PValue, "fdr"), comparison_name
	)][, DE := case_when(
		log2FC >= log2fc_cutoff & FDR <= fdr_cutoff ~ "up",
		log2FC <= -log2fc_cutoff & FDR <= fdr_cutoff ~ "down",
		TRUE ~ "unchanged"
	)]

	## Merge in the annotation information from the original matrix.
	setkey(diff_expression, "FID")
	setkey(anno_data, "FID")
	comparison_name <- str_c(compare_groups[1], "_vs_", compare_groups[2])
	diff_expression <- diff_expression[anno_data, nomatch = 0]

	## Add differential expression data back to object.
	if (data_type == "tss") {
		experiment@diff_features$TSSs[[comparison_name]] <- diff_expression
	} else if (data_type == "tsr") {
		experiment@diff_features$TSRs[[comparison_name]] <- diff_expression
	} else if (data_type == "tss_features") {
		experiment@diff_features$TSS_features[[comparison_name]] <- diff_expression
	} else if (data_type == "tsr_features") {
		experiment@diff_features$TSR_features[[comparison_name]] <- diff_expression
	}

	return(experiment)
}

#' DE MA Plot
#'
#' Generate MA plot for differential TSRs or Genes (RNA-seq)
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr case_when mutate
#'
#' @param experiment tsrexplorer object
#' @param de_comparisons Which differential expression comparisons to plot
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param ncol Number of columns for the facets
#' @param ... Arguments passed to geom_point
#'
#' @return ggplot2 object of differential TSRs volcano plot.
#'
#' @rdname plot_volcano-function
#'
#' @export

plot_ma <- function(
	differential_expression,
	data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	de_comparisons = "all", ncol = 1, ...
){

	## Get differential expression tables.
	if (data_type == "tss") {
		de_samples <- experiment@diff_features$TSSs
	} else if (data_type == "tsr") {
		de_samples <- experiment@diff_features$TSRs
	} else if (data_type == "tss_features") {
		de_samples <- experiment@diff_features$TSS_features
	} else if (data_type == "tsr_features") {
		de_samples <- experiment@diff_features$TSR_features
	}

	if (de_comparisons == "all") {
		de_samples <- discard(de_samples, names(de_samples) %in% c("model", "design"))
	} else {
		de_samples <- de_samples[de_comparisons]
	}

	de_samples <- bind_rows(de_samples)
	de_samples <- de_samples[, .(sample, FID, log2FC, logCPM, DE)]
	de_samples[, DE := factor(DE, levels = c("up", "unchanged", "down"))]

	## MA plot of differential expression
	p <- ggplot(de_samples, aes(x = logCPM, y = log2FC, color = DE)) +
		geom_point(...) +
		theme_bw() +
		scale_color_viridis_d() +
		facet_wrap(~ sample, ncol = ncol, scales = "free")

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
