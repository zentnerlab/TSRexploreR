
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
#' @param experiment tsrexplorer object with TMM-normalized counts
#' @param data_type Whether TSSs, TSRs, or feature counts should be analyzed
#' @param samples Vector of sample names to analyze
#' @param groups Vector of groups in correct factor notation
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

	## Set sample design.
	sample_design <- model.matrix(~ 0 + design[["groups"]])

	## Create DE model.
	fitted_model <- assay(sample_data, "counts") %>%
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

#' AnalyzeDifferential Expression
#'
#' Find differential TSSs, TSRs, or features from edgeR model
#'
#' @import tibble
#' @importFrom edgeR glmQLFTest
#' @importFrom dplyr pull mutate
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
#'
#' @param experiment tsrexplorer object with edgeR differential expression model from fit_edger_model
#' @param data_type Whether the input was generated from TSSs, TSRs, or features
#' @param compare_groups Vector of length two of the two groups from which to find differential TSRs
#' @param fdr_cutoff FDR cutoff
#' @param log2fc_cutoff Log2 fold change cutoff
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
		original_ranges <- rowRanges(experiment@counts$TSSs$matrix)
	} else if (data_type == "tsr") {
		edger_model <- experiment@diff_features$TSRs$model
		edger_design <- experiment@diff_features$TSRs$design
		original_ranges <- rowRanges(experiment@counts$TSRs$matrix)
	} else if (data_type == "tss_features") {
		edger_model <- experiment@diff_features$TSS_features$model
		edger_design <- experiment@diff_features$TSS_features$design
	} else if (data_type == "tsr_features") {
		edger_model <- experiment@diff_features$TSR_features$model
		edger_design <- experiment@diff_features$TSR_features$design
	}

	## Set up contrasts.
	comparison_levels <- edger_design %>%
		pull(groups) %>%
		levels

	comparison_contrast <- comparison_levels %>%
		map_dbl(function(x) {
			if (x == compare_groups[1]) {
				return_value <- -1
			} else if (x == compare_groups[2]) {
				return_value <- 1
			} else {
				return_value <- 0
			}
			return(return_value)
		})

	## Differential expression analysis.
	diff_expression <- glmQLFTest(edger_model, contrast = comparison_contrast)
	diff_expression <- as.data.table(diff_expression$table, keep.rownames = "FHASH")

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
	original <- as.data.table(original_ranges)
	original[, FHASH := names(original_ranges)]

	comparison_name <- str_c(compare_groups[1], "_vs_", compare_groups[2])
	diff_expression <- merge(diff_expression, original, by = "FHASH")

	diff_expression <- sort(makeGRangesFromDataFrame(diff_expression, keep.extra.columns = TRUE))
	diff_expression <- as.data.table(diff_expression)
	

	## Add differential expression data back to tsrexplorer object.
	if (data_type == "tss") {
		experiment@diff_features$TSSs$results[[comparison_name]] <- diff_expression
	} else if (data_type == "tsr") {
		experiment@diff_features$TSRs$results[[comparison_name]] <- diff_expression
	} else if (data_type == "tss_features") {
		experiment@diff_features$TSS_features$results[[comparison_name]] <- diff_expression
	} else if (data_type == "tsr_features") {
		experiment@diff_features$TSR_features[[comparison_name]] <- diff_expression
	}

	return(experiment)
}

#' DE Table
#'
#' Output a table with differential features
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The name of the DE comparison
#' @param de_type A single value or combination of 'up, 'unchanged', and/or 'down' (qq a list?)
#'
#' @rdname de_table-function
#' @export

de_table <- function(
	experiment, data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	de_comparisons = "all", de_type = c("up", "unchanged", "down")
) {
	## Grab tables.
	de_tables <- experiment %>%
		extract_de(data_type, de_comparisons) %>%
		bind_rows

	## Filter tables.
	de_tables <- de_tables[DE %in% de_type]

	## Return tables.
	return(de_tables)
}
