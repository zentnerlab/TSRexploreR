
#' Genes/Transcripts Detected
#'
#' Get the number of genes or transcripts detected
#'
#' @import tibble
#' @importFrom dplyr bind_rows full_join rename filter select mutate distinct case_when count
#' @importFrom purrr map
#' @importFrom magrittr extract %>%
#' @importFrom SummarizedExperiment assay rowRanges SummarizedExperiment
#' @importFrom S4Vectors DataFrame "metadata<-"
#'
#' @param experiment tsrexplorer object with annotated TSSs or TSRs
#' @param samples Either 'all' or vector of sample names
#' @param data_type Whether TSSs or TSRs should be analyzed
#' @param threshold The number of reads required in a TSS or TSR to avoid filtering
#'
#' @return tibble of detected feature numbers
#'
#' @rdname detect_features-function
#'
#' @export

detect_features <- function(
	experiment,
	samples = "all",
	data_type = c("tss", "tsr"),
	threshold = 1,
	dominant = FALSE
) {
	
	## Get sample data.
	sample_data <- experiment %>%
		extract_counts(data_type, samples) %>%
		bind_rows(.id = "sample")

	setnames(
		sample_data, old = ifelse(
			experiment@settings$annotation[, feature_type] == "transcript",
			"transcriptId", "geneId"
		),
		new = "feature"
	)

	## prepare sample data.
	keep_cols <- c("sample", "feature", "simple_annotations")
	if (dominant) keep_cols <- c(keep_cols, "dominant")

	sample_data <- sample_data[
		score >= threshold,
		..keep_cols
	]

	## Get feature counts.
	sample_data <- sample_data[,
		.(promoter = any(simple_annotations == "Promoter")),
		by = .(sample, feature)
	][,
		.(with_promoter_tss = sum(promoter), no_promoter_tss = .N - sum(promoter), total = .N),
		by = sample
	]

	## Create DataFrame to export.
	detected_features <- DataFrame(sample_data)

	metadata(detected_features)$threshold <- threshold
	metadata(detected_features)$data_type <- data_type
	metadata(detected_features)$dominant <- dominant
	metadata(detected_features)$feature_type <- experiment@settings$annotation[, feature_type]

	return(detected_features)
}

#' Plot Detected Features
#'
#' Plot number of features detected per sample.
#'
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom stringr str_c str_to_title
#'
#' @param detected_features Tibble of detected feature counts from detect_features
#' @param ... Arguments passed to geom_col
#'
#' @return ggplot2 object of detected feature counts
#'
#' @rdname plot_detected_features-function
#'
#' @export

plot_detected_features <- function(detected_features, ...) {
	
	## Get some info from DataFrame.
	feature_type <- str_to_title(metadata(detected_features)$feature_type)	

	## Prepare data for plotting.
	plot_data <- as.data.table(detected_features)
	plot_data[, total := NULL]	
	plot_data <- melt(plot_data,
		measure.vars = c("with_promoter_tss", "no_promoter_tss"),
		variable.name = "count_type",
		value.name = "feature_count"
	)

	## Plot data.
	p <- ggplot(plot_data, aes(x = sample, y = feature_count, fill = fct_rev(factor(count_type)))) +
		geom_col(position = "stack", ...) +
		theme_bw() +
		scale_fill_viridis_d(name = "Feature Type", direction = -1) +
		ylim(c(0, NA)) +
		ylab(str_c(feature_type, "Count", sep = " ")) +
		xlab("Sample") +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1)
		)
		
}
