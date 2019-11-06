
#' Genes/Transcripts Detected
#'
#' Get the number of genes or transcripts detected
#'
#' @import tibble
#' @importFrom dplyr rename bind_rows select mutate case_when between distinct full_join count
#' @importFrom purrr map
#'
#' @param experiment tsrexplorer object with annotated TSSs or TSRs
#' @param samples Either 'all' or vector of sample names
#' @param data_type Whether TSSs or TSRs should be analyzed
#' @param feature_type Whether the TSSs and TSRs were annotated relative to genes or transcripts
#' @param upstream Number of bases upstream of annotated TSS to consider as part of promoter
#' @param downstream Number of bases downstream of annotated TSS to consider as part of the promoter
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
	feature_type = c("gene", "transcript"),
	upstream = 1000,
	downstream = 1000,
	threshold = 1
) {
	
	## Ensure appropriate sample names if "all" selected.
	if (data_type == "tss") {
		if (samples == "all") samples <- names(experiment@annotated$TSSs)
		sample_data <- experiment@annotated$TSSs[samples]
	} else if (data_type == "tsr") {
		if (samples == "all") samples <- names(experiment@annotated$TSRs)
		sample_data <- experiment@annotated$TSRs[samples] %>%
			map(~ rename(., score = nTAGs))
	}

	## Pull and combine chosen sample data.
	sample_data <- bind_rows(sample_data, .id = "sample")

	## Pick feature type to use.
	if (feature_type == "gene") {
		sample_data <- rename(sample_data, feature = geneId)
	} else if (feature_type == "transcript") {
		sample_data <- rename(sample_data, feature = transcriptId)
	}

	## Get feature counts.
	sample_data <- sample_data %>%
		filter(score >= threshold) %>%
		select(sample, feature, distanceToTSS) %>%
		mutate(annotation = case_when(
			between(distanceToTSS, -upstream, downstream) ~ "promoter",
			TRUE ~ "other"
		)) %>%
		select(-distanceToTSS)

	total_features <- sample_data %>%
		distinct(sample, feature) %>%
		count(sample, name = "total_features")

	promoter_proximal_features <- sample_data %>%
		filter(annotation == "promoter") %>%
		distinct(sample, feature) %>%
		count(sample, name = "promoter_proximal_features")

	detected_features <- full_join(total_features, promoter_proximal_features, by = "sample")

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
#'
#' @param detected_features Tibble of detected feature counts from detect_features
#' @param ncol Number of columens to plot data to
#'
#' @return ggplot2 object of detected feature counts
#'
#' @rdname plot_detected_features-function
#'
#' @export

plot_detected_features <- function(detected_features, ncol = 1) {
	
	## Prepare data for plotting.
	plot_data <- detected_features %>%
		gather(key = "feature_type", value = "feature_number", -sample) %>%
		mutate(feature_type = factor(feature_type, levels = c("total_features", "promoter_proximal_features")))

	## Plot data.
	p <- ggplot(plot_data, aes(x = feature_type, y = feature_number, fill = feature_type)) +
		geom_col() +
		theme_bw() +
		scale_fill_viridis_d(name = "Feature Type") +
		facet_wrap(~ sample, ncol = ncol) +
		ylim(c(0, NA)) +
		ylab("Feature Count") +
		theme(
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			axis.title.x = element_blank()
		)
		
}
