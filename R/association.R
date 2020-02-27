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

#' Associate TSSs
#'
#' Associate TSSs with TSRs
#'
#' @importFrom purrr iwalk
#'
#' @param experiment tsrexplorer object
#' @param use_sample_sheet Whether to associate TSSs with their matching TSR in the sample sheet
#' @param sample_list If 'use_sample_sheet' is FALSE, provide list with TSR as name and vector of TSSs to associate with it.
#'
#' @rdname associate_with_tsr-function
#' @export

associate_with_tsr <- function(experiment, use_sample_sheet = TRUE, sample_list = NA) {

	## Grab sample sheet if being used.
	if (use_sample_sheet) {
		samples <- experiment@meta_data$sample_sheet
		sample_list <- as.list(samples$tss_name)
		sample_list <- set_names(sample_list, samples$tsr_name)
	}
	
	## Associate TSSs with TSRs.
	associated_TSSs <- imap(sample_list, function(tss_names, tsr_name) {

		# Make GRanges of TSSs.
		tss_set <- extract_counts(experiment, "tss", tss_names) %>%
			bind_rows(.id = "sample")
		
		tss_gr <- makeGRangesFromDataFrame(tss_set, keep.extra.columns = TRUE)

		# Make GRanges of TSRs.
		tsr_set <- extract_counts(experiment, "tsr", tsr_name) %>%
			bind_rows(.id = "tsr_sample")
		
		tsr_gr <- tsr_set
		setnames(tsr_gr, old = "FID", new = "TSR_FID")
		tsr_gr <- tsr_gr[, .(seqnames, start, end, strand, TSR_FID)]
		tsr_gr <- makeGRangesFromDataFrame(tsr_gr, keep.extra.columns = TRUE)

		# Annotate TSS to overlapping TSR.
		overlapping <- as.data.table(findOverlapPairs(tss_gr, tsr_gr))
		setnames(
			overlapping,
			old = c(
				"first.X.seqnames", "first.X.start", "first.X.end",
				"first.X.strand", "second.TSR_FID", "first.sample"
			),
			new = c("seqnames", "start", "end", "strand", "TSR_FID", "sample")
		)
		overlapping <- overlapping[, .(seqnames, start, end, strand, TSR_FID, sample)]

		# Prepare tsr_set to merge into overalpping tss.
		tsr_clean <- copy(tsr_set)
		tsr_clean[, tsr_coords := str_c(seqnames, start, end, strand, sep = ":")]
		setnames(
			tsr_clean, old = c("score", "width"), new = c("tsr_score", "tsr_width")
		)
		tsr_clean[, c("seqnames", "start", "end", "strand") := NULL]

		# Merge overlapping back into TSS.
		key_ids <- c("sample", "seqnames", "start", "end", "strand")
		setkeyv(overlapping, key_ids)
		setkeyv(tss_set, key_ids)
		tss_set <- merge(tss_set, overlapping, all.x = TRUE)

		# Merge TSRs into TSSs.
		setkey(tsr_clean, "TSR_FID")
		setkey(tss_set, "TSR_FID")
		merged_set <- merge(tss_set, tsr_clean, all.x = TRUE)

		merged_set <- merged_set %>%
			split(.$sample) %>%
			map(function(x) {
				x <- x[order(FID)]
				x[, sample := NULL]
				x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)			
				x <- as.data.table(x)
				return(x)
			})
	})

	## Add TSSs back to tsrexplorer object.
	associated_TSSs <- purrr::reduce(associated_TSSs, c)
	
	experiment@counts$TSSs$raw <- associated_TSSs
	return(experiment)
}
