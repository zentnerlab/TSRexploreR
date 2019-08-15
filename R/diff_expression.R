
#' Build edgeR Model of Differential TSRs
#'
#' Find differential TSRs
#'
#' @import tibble
#' @importFrom edgeR DGEList filterByExpr calcNormFactors cpm estimateDisp glmQLFit
#' @importFrom dplyr select_at
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer objet after TMM normalization.
#' @param samples Vector of sample names to analyze.
#' @param groups Vector of groups in correct factor notation.
#'
#' @return DGEList object of differentially epxressed genes.with fitted DE model.
#'
#' @export
#' @rdname fit_edger_model-function

fit_edger_model <- function(experiment, samples = c(), groups = c()) {

	## Select samples and turn to count matrix.
	selected_samples <- experiment@raw_counts$TSRs %>%
		select_at(c("TSR_name", samples)) %>%
		column_to_rownames("TSR_name") %>%
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

#' Find Differential TSRs
#'
#' Find differential TSRs from edgeR model.
#'
#' @import tibble
#' @importFrom edgeR glmQLFTest
#' @importFrom dplyr pull
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#'
#' @param fit_edger_model edgeR differential expression model from fit_edger_model
#' @param comparison Vector of length two of the two groups to find differential TSRs
#'
#' @return tibble of differential TSRs
#'
#' @export
#' @rdname differential_tsrs-function

differential_tsrs <- function(fit_edger_model, comparisons = c()) {
	
	## Set up contrasts.
	comparison_contrast <- fit_edger_model$samples %>%
		pull(group) %>%
		levels %>%
		as.numeric %>%
		length %>%
		numeric(length = .)

	comparison_contrast[comparisons[1]] <- 1
	comparison_contrast[comparisons[2]] <- -1

	## Differential expression
	diff_expression <- glmQLFTest(fit_edger_model, contrast = comparison_contrast)

	## Prepare tibble for export
	diff_expression <- diff_expression$table %>%
		as_tibble(.name_repair = "universal", rownames = "position") %>%
		separate(position, into = c("chr", "start", "end", "strand"), sep = "_") %>%
		mutate(FDR = p.adjust(PValue, method = "fdr"))

	return(diff_expression)
}

#' Annotate Differential TSRs
#'
#' Annotate Differential TSRs to nearest gene or transcript.
#'
#' @import tibble
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom magrittr %>%
#'
#' @param differential_tsrs Tibble of differential TSRs from differential_tsrs
#' @param annotation_file GTF genomic annotation file
#' @param feature_type Whether to annotate TSRs relative to genes to transcripts
#' @param upstream Bases upstream of TSS
#' @param downstream Bases downstream of TSS
#'
#' @return Tibble of annotated differential TSRs
#'
#' @export
#' @rdname annotate_differential_tsrs-function

annotate_differential_tsrs <- function(
	differential_tsrs, annotation_file,
	feature_type = c("gene", "transcript"),
	upstream = 1000, downstream = 100
) {
	## Load genome annotation file as TxDb.
	genome_annotation <- makeTxDbFromGFF(annotation_file, "gtf")

	## Annotate differential TSRs.
	annotated_diff_tsrs <- differential_tsrs %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
		annotatePeak(
			tssRegion = c(-upstream, downstream),
			TxDb = genome_annotation,
			sameStrand = TRUE,
			level = feature_type
		) %>%
		as_tibble(.name_repair = "universal")

	return(annotated_diff_tsrs)
}
