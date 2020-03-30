
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
#' @param condition_data Apply conditions to data (supports filtering and quantiles/grouping)
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
	threshold = NA,
	dominant = FALSE,
	condition_data = NA
) {
	
	## Get sample data.
	sample_data <- extract_counts(experiment, data_type, samples)

	## Initial sample processing.
	if (dominant | !is.na(threshold)) {
		sample_data <- preliminary_filter(sample_data, dominant, threshold)
	}

	## Apply data conditioning if set.
	if (!is.na(condition_data)) {
		sample_data <- do.call(group_data, c(list(signal_data = sample_data), condition_data))
	}

	## Rename feature column.
	walk(sample_data, function(x) {
		setnames(
			x, old = ifelse(
				experiment@settings$annotation[, feature_type] == "transcript",
				"transcriptId", "geneId"
			),
			new = "feature"
		)
	})

	## Get feature counts.
	groupings <- any(names(data_group) %in% c("quantile_by", "grouping"))
	sample_data <- rbindlist(sample_data, idcol = "sample")
	
	if (groupings) {
		sample_data <- sample_data[,
			.(grouping, promoter = any(simple_annotations == "Promoter")),
			by = .(sample, feature)
		][,
			.(with_promoter = sum(promoter), without_promoter = .N - sum(promoter), total = .N),
			by = .(sample, grouping)
		]
	} else {
		sample_data <- sample_data[,
			.(promoter = any(simple_annotations == "Promoter")),
			by = .(sample, feature)
		][,
			.(with_promoter = sum(promoter), without_promoter = .N - sum(promoter), total = .N),
			by = sample
		]
	}

	## Create DataFrame to export.
	detected_features <- DataFrame(sample_data)

	metadata(detected_features)$threshold <- threshold
	metadata(detected_features)$data_type <- data_type
	metadata(detected_features)$dominant <- dominant
	metadata(detected_features)$groupings <- groupings
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
	grouping_status <- metadata(detected_features)$groupings

	## Prepare data for plotting.
	plot_data <- as.data.table(detected_features)
	plot_data[, total := NULL]	
	plot_data <- melt(plot_data,
		measure.vars = c("with_promoter", "without_promoter"),
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

	if (grouping_status) {
		p <- p + facet_grid(fct_rev(factor(grouping)) ~ .)
	}
		
	return(p)
}
