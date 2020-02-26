
#' Format Counts
#'
#' Format counts for TSSs, TSRs, and/or features.
#'
#' @importFrom magrittr set_colnames set_rownames
#'
#' @param experiment tsrexplorer object
#' @param data_type 'tss' or 'tsr'
#' @param samples Sample names or merged replicate IDs to process
#'
#' @rdname format_counts-function
#' @export

format_counts <- function(experiment, data_type = c("tss", "tsr"), samples = "all") {

	## Grab appropriate samples and generate raw count matrices.
	if (data_type == "tss") {

		## Grab selected_samples.
		select_samples <- tss_experiment(experiment)
		if (samples != "all") select_samples <- select_samples[samples]

	} else if (data_type == "tsr") {

		## Grab selected samples.
		select_samples <- tsr_experiment(experiment)
		if (samples != "all") select_samples <- select_samples[samples]

	}

	## Turn counts into data.table
	if (data_type %in% c("tss", "tsr")) {
		raw_counts <- map(select_samples, function(x) {
			x <- as.data.table(x)
			x[, FID := seq_len(nrow(x))]
			return(x)
		})
	}

	## Place counts in proper object slot.
	if (data_type == "tss") {
		experiment@counts$TSSs$raw <- raw_counts
	} else if (data_type == "tsr") {
		experiment@counts$TSRs$raw <- raw_counts
	}

	return(experiment)
}

#' Feature Counts
#'
#' Add feature counts
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss' or 'tsr'
#'
#' @rdname count_features-function
#' @export

count_features <- function(experiment, data_type = c("tss", "tsr")) {
	
	## Get information on whether annotation was by gene or transcript.
	anno_type <- ifelse(
		experiment@settings$annotation[, feature_type] == "transcript",
		"transcriptId", "geneId"
	)

	## Extract appropriate counts.
	sample_data <- extract_counts(experiment, data_type, "all")

	## Get feature counts.
	sample_data <- sample_data %>%
		map(function(x) {
			feature_scores <- x[
				!simple_annotations %in% c("Downstream", "Intergenic"),
				.(feature_score = sum(score)),
				by = eval(anno_type)
			]

			setkeyv(x, anno_type)
			setkeyv(feature_scores, anno_type)
			x <- merge(x, feature_scores, all.x = TRUE)
			setcolorder(x, c(discard(colnames(x), ~ . == anno_type), anno_type))
			
			return(x)
		})

	## Get feature counts per sample.
	counts <- sample_data %>%
                map(function(x) {
                        setnames(x, old = anno_type, new = "feature")
                        x <- unique(x[, .(feature, feature_score)])
                        x[, feature_score := ifelse(is.na(feature_score), 0, feature_score)]
                        return(x)
                })

	## Store counts in appropriate slots.
	walk(counts, ~ setnames(., old = c("feature", "feature_score"), new = c(anno_type, "score")))
	walk(sample_data, ~ setnames(., old = "feature", new = anno_type))

	if (data_type == "tss") {
		experiment@counts$TSSs$raw <- sample_data
		experiment@counts$TSS_features$raw <- counts
	} else {
		experiment@counts$TSRs$raw <- sample_data
		experiment@counts$TSR_features$raw <- exp_counts
	}

	return(experiment)
}

#' Count Matrix
#'
#' Generate count matrices
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param samples The samples to turn into a count matrix
#'
#' @rdname count_matrix-function
#' @export

count_matrix <- function(
	experiment, data_type = c("tss", "tsr", "tss_features", "tsr_features"),
	samples = "all"
) {
	## Extract counts.
	select_samples <- extract_counts(experiment, data_type, samples)

	if (data_type %in% c("tss", "tsr")) {
		select_samples <- map(select_samples, function(x) {
			x <- x[, .(seqnames, start, end, strand, score)]
			return(x)
		})
	}

	## Turn into matrix.
	if (data_type == "tss") {

		## Merge overlapping TSSs
		select_samples <- bind_rows(select_samples, .id = "sample")
		select_samples <- dcast(
			select_samples,
			seqnames + start + end + strand ~ sample,
			fill = 0
		)

	} else if (data_type == "tsr") {

                ## Merge overlapping TSRs to get consensus.
                tsr_consensus <- select_samples %>%
                        map(makeGRangesFromDataFrame) %>%
			purrr::reduce(c) %>%
                        GenomicRanges::reduce(ignore.strand = FALSE)

                ## Create raw count matrix.
                raw_matrix <- select_samples %>%
                        map(
                                ~ makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>%
				findOverlapPairs(query = tsr_consensus, subject = .) %>%
                                as.data.table
                        ) %>%
                        bind_rows(.id = "sample")

                setnames(
                        raw_matrix,
                        old = c(
                                "first.seqnames", "first.start", "first.end",
                                "first.strand", "second.X.score"
                        ),
                        new = c("seqnames", "start", "end", "strand", "score")
                )

                raw_matrix <- raw_matrix[,
                        .(score = sum(score)),
                        by = .(sample, seqnames, start, end, strand)
                ]

                select_samples <- dcast(raw_matrix, seqnames + start + end + strand ~ sample, fill = 0)

	} else if (data_type %in% c("tss_features", "tsr_features")) {
		## Get annotation type.
		anno_type <- experiment@settings$annotation[["feature_type"]]
		anno_type <- ifelse(anno_type == "gene", "geneId", "transcriptId")

		## Change feature counts to matrix
		select_samples <- bind_rows(select_samples, .id = "sample")
		select_samples <- dcast(
			select_samples,
			as.formula(str_c(anno_type, " ~ sample")),
			fill = 0
		)

	}

	## Create SummarizedExperiments.
	if (data_type %in% c("tss", "tsr")) {
		## Prepare data for RangedSummarizedExperiment
		row_ranges <- makeGRangesFromDataFrame(select_samples)

		select_samples[, c("seqnames", "start", "end", "strand") := NULL]
		select_samples <- as.matrix(select_samples)

        	col_data <- DataFrame(
        		samples = colnames(select_samples),
                	row.names = colnames(select_samples)
        	)
	
		## Make RangedSummarizedExperiment.
		rse <- SummarizedExperiment(
			assays = list(counts = select_samples),
			rowRanges = row_ranges,
			colData = col_data
		)
	} else if (data_type %in% c("tss_features", "tsr_features")) {
		## Prepare data for SummarizedExperiment.
		row_data <- select_samples[, 1]
		
		select_samples[, eval(anno_type) := NULL]
		select_samples <- as.matrix(select_samples)

		col_data <- DataFrame(
			samples = colnames(select_samples),
			row.names = colnames(select_samples)
		)

		## Make SummarizedExperiment.
		rse <- SummarizedExperiment(
			assays = list(counts = select_samples),
			rowData = row_data,
			colData = col_data
		)
	}

	## Add RangedSummarizedExperiment back to trexplorer object.
	if (data_type == "tss") {
		experiment@counts$TSSs$matrix <- rse
	} else if (data_type == "tsr") {
		experiment@counts$TSRs$matrix <- rse
	} else if (data_type == "tss_features") {
		experiment@counts$TSS_features$matrix <- rse
	} else if (data_type == "tsr_features") {
		experiment@counts$TSR_features$matrix <- rse
	}

	return(experiment)	
}





