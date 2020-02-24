
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

		## Create raw count matrix.
		raw_matrix <- select_samples %>%
			map(as.data.table) %>%
			bind_rows(.id = "sample")		
		raw_matrix <- dcast(raw_matrix, seqnames + start + end + strand ~ sample, fill = 0)

	} else if (data_type == "tsr") {

		## Grab selected samples.
		select_samples <- tsr_experiment(experiment)
		if (samples != "all") select_samples <- select_samples[samples]

		## Merge overlapping TSRs to get consensus.
		tsr_consensus <- select_samples %>%
			purrr::reduce(c) %>%
			GenomicRanges::reduce(ignore.strand = FALSE)

		## Create raw count matrix.
		raw_matrix <- select_samples %>%
			map(
				~ findOverlapPairs(query = tsr_consensus, subject = .) %>%
				as.data.table
			) %>%
			bind_rows(.id = "sample")
		raw_matrix <- setnames(
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
		raw_matrix <- dcast(raw_matrix, seqnames + start + end + strand ~ sample, fill = 0)
	
	}

	## Create RangedSummarizedExperiment for regular counts.
	if (data_type %in% c("tss", "tsr")) {
		raw_counts <- select_samples %>%
			imap(function(gr, sample_name) {
				count_data <- gr %>%
					score(.) %>%
					as.matrix %>%
					set_colnames(sample_name)

				row_data <- gr
				score(row_data) <- NULL
				row_data$FID <- sprintf("FID%s", seq_len(length(row_data)))
				col_data <- DataFrame(sample = sample_name)

				raw_exp <- SummarizedExperiment(
					assays = list(raw = count_data),
					rowRanges = row_data,
					colData = col_data
				)
				return(raw_exp)
			})
	}

	## Create RangedSummarizedExperiment for count matrices.
	if (data_type %in% c("tss", "tsr")) {
		count_matrix <- as.matrix(raw_matrix[,-1:-4])

		row_data <- makeGRangesFromDataFrame(raw_matrix)
		row_data$FID <- sprintf("FID%s", seq_len(length(row_data)))
	}
	
	col_data <- DataFrame("sample" = colnames(count_matrix), row.names = colnames(count_matrix))
	
	matrix_counts <- SummarizedExperiment(
		assays = list(counts = count_matrix),
		rowRanges = row_data,
		colData = col_data
	)

	## Place counts in proper object slot.
	if (data_type == "tss") {
		experiment@counts$TSSs <- list(
			"raw" = raw_counts,
			"matrix" = matrix_counts
		)
	} else if (data_type == "tsr") {
		experiment@counts$TSRs <- list(
			"raw" = raw_counts,
			"matrix" = matrix_counts
		)
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
	sample_data <- experiment %>%
		extract_counts(data_type, "all") %>%
		map(as.data.table)

	## Get feature counts.
	sample_data <- sample_data %>%
		map(function(x) {
			x <- x[
				!simple_annotations %in% c("Downstream", "Intergenic"),
				.(score = sum(score)),
				by = eval(anno_type)
			]
			setnames(x, old = anno_type, new = "feature")
			return(x)
		})

	## Turn feature counts into SummarizedExperiments.
	exp_counts <- sample_data %>%
		imap(function(x, y) {
			setnames(x, old = "score", new = y)
			mat <- x %>%
				column_to_rownames("feature") %>%
				as.matrix

			row_data <- DataFrame("feature" = rownames(mat))
			col_data <- DataFrame("sample" = colnames(mat))

			se <- SummarizedExperiment(
				assay = list("raw" = mat),
				colData = col_data,
				rowData = row_data
			)
			return(se)
		})

	## Create count matrix.
	count_matrix <- sample_data %>%
		imap(~ setnames(.x, old = .y, new = "score")) %>%
		bind_rows(.id = "sample") %>%
		dcast(feature ~ sample, fill = 0) %>%
		column_to_rownames("feature") %>%
		as.matrix

	## Turn count matrix into SummarizedExperiment.
	row_data <- DataFrame(feature = rownames(count_matrix))
	col_data <- DataFrame(sample = colnames(count_matrix))

	exp_matrix <- SummarizedExperiment(
		assay= list("counts" = count_matrix),
		rowData = row_data,
		colData = col_data
	)

	## Store matrix and counts in appropriate slots.
	if (data_type == "tss") {
		experiment@counts$TSS_features$raw <- exp_counts
		experiment@counts$TSS_features$matrix <- exp_matrix
	} else {
		experiment@counts$TSR_features$raw <- exp_counts
		experiment@counts$TSR_features$matrix <- exp_matrix
	}

	return(experiment)
}

#' Merge Samples
#'
#' Merge replicates or selected samples.
#'
#' @param experiment tsrexplorer object
#' @param data_type Either 'tss' or 'tsr'
#' @param merge_replicates If 'TRUE' replicate groups will be merged
#' @param sample_list If merge_replicates is set to 'FALSE',
#' specify what samples to merge in list format.
#'
#' @rdname merge_samples-function
#' @export

merge_samples <- function(
	experiment, data_type = c("tss", "tsr"),
	merge_replicates = TRUE, sample_list = NA
) {
	
	## Prepare what samples will be merged.
	if (merge_replicates) {
		samples <- experiment@meta_data$sample_sheet[type == data_type] %>%
			split(.$replicate_id) %>%
			map(~ pull(., name))
	}

	## Get feature sets to be merged.
	if (data_type == "tss") {
		selected_samples <- experiment@experiment$TSSs[unlist(samples)]
	} else if (data_type == "tsr") {
		selected_samples <- experiment@experiment$TSRs[unlist(samples)]
	}

	selected_samples <- map(selected_samples, as.data.table)

	## Merge feature sets.
	if (data_type == "tss") {
		merged_samples <- samples %>%
			map(function(sample_group) {
				merged <- bind_rows(selected_samples[sample_group])
				merged <- as.data.table(merged)
				merged[, score := sum(score), by = .(seqnames, start, end, strand)]
				return(merged)
			})
	} else if (data_type == "tsr") {
		merged_samples <- samples %>%
			map(function(sample_group) {
				merged <- selected_samples[sample_group] %>%
					map(~ makeGRangesFromDataFrame(., keep.extra.columns = TRUE))
				
				tsr_consensus <- merged %>%
					purrr::reduce(c) %>%
					GenomicRanges::reduce(ignore.strand = FALSE)

				merged <- merged %>%
					map(
						~ findOverlapPairs(query = tsr_consensus, subject = .) %>%
						as.data.table
					) %>%
					bind_rows

				setnames(
					merged,
					old = c(
						"first.seqnames", "first.start", "first.end",
						"first.strand", "second.X.score"
					),
					new = c("seqnames", "start", "end", "strand", "score")
				)

				merged <- merged[,
					.(score = sum(score)),
					by = .(seqnames, start, end, strand)
				]
			})
	}

	## Convert merged samples to GRanges.
	merged_samples <- map(merged_samples, ~ makeGRangesFromDataFrame(., keep.extra.columns = TRUE))

	## Return merged samples.
	if (data_type == "tss") {
		experiment@experiment$TSSs <- c(experiment@experiment$TSSs, merged_samples)
	} else if (data_type == "tsr") {
		experiment@experiment$TSRs <- c(experiment@experiment$TSRs, merged_samples)
	}

	return(experiment)
}
