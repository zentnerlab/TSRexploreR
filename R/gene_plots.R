
#' Genes/Transcripts Detected
#'
#' Get the number of genes or transcripts detected
#'
#' @importFrom purrr walk
#'
#' @param experiment tsrexplorer object with annotated TSSs or TSRs
#' @param samples Either 'all' or vector of sample names
#' @param data_type Whether TSSs or TSRs should be analyzed
#' @param threshold The number of reads required in a TSS or TSR to avoid filtering
#' @param dominant Whether to consider only the dominant TSS or TSR
#' @param condition_data Apply conditions to data (supports filtering and quantiles/grouping)
#'
#' @return DataFrame of detected feature numbers.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data = annotation,
#'   data_type = "tss", feature_type = "transcript"
#' )
#' detected <- detect_features(tsre_exp, data_type = "tss")
#'
#' @seealso
#' \code{\link{annotate_features}} to annotate the TSSs and TSRs.
#' \code{\link{plot_detected_features}} to plot detected features.
#'
#' @rdname detect_features-function
#' @export

detect_features <- function(
	experiment,
	samples = "all",
	data_type = c("tss", "tsr"),
	threshold = NA,
	dominant = FALSE,
	condition_data = NA
) {

	## Check inputs.
	if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")

	if (!is(samples, "character")) stop("samples must be a character vecotor")

        if (!is(data_type, "character")) stop("data_type must be a character")
        if (length(data_type) > 1) stop("data_type must be a character")
        data_type <- str_to_lower(data_type)
        if (!data_type %in% c("tss", "tsr")) stop("data_type must be 'tss' or 'tsr'")

        if (!is.na(threshold) & !is(threshold, "numeric")) stop("threshold must be a positive integer")
        if (!is.na(threshold) & threshold %% 1 != 0) stop("threshold must be a positive integer")
        if (!is.na(threshold) & threshold < 1) stop("threshold must be greater than or equal to 1")

	if (!is(dominant, "logical")) stop("dominant must be logical")

	if (!is.na(condition_data) & !is(condition_data, "list")) {
		stop("condition_data must be a list of values")
	}
	
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
	groupings <- any(names(condition_data) %in% c("quantile_by", "grouping"))
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
#' @importFrom stringr str_to_title
#'
#' @param detected_features DataFrame of detected feature counts from detect_features
#' @param ... Arguments passed to geom_col
#'
#' @return ggplot2 object of detected feature counts
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data = annotation,
#'   data_type = "tss", feature_type = "transcript"
#' )
#' detected <- detect_features(tsre_exp, data_type = "tss")
#' plot_detected_features(detected)
#'
#' @seealso
#' \code{\link{annotate_features}} to annotate the TSSs or TSRs.
#' \code{\link{detect_features}} to first detect feature numbers.
#'
#' @rdname plot_detected_features-function
#' @export

plot_detected_features <- function(
	detected_features,
	...
) {

	## Check inputs.
	if (!is(detected_features, "DataFrame")) stop("detected_features must be a DataFrame")
	
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
