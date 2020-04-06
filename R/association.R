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
	merge_replicates = FALSE, sample_list = NA
) {
	
	## Prepare what samples will be merged.
	if (merge_replicates) {
		data_column <- ifelse(data_type == "tss", "tss_name", "tsr_name")

		sample_list <- experiment@meta_data$sample_sheet[type == data_type] %>%
			split(.$replicate_id) %>%
			map(~ pull(., name))
	}

	## Merge feature sets.
	if (data_type == "tss") {
		merged_samples <- map(sample_list, function(sample_group) {
			select_samples <- experiment@experiment$TSSs[sample_group]
			select_samples <- map(select_samples, as.data.table)

			merged <- rbindlist(select_samples)
			merged[, .(score = sum(score)), by = .(seqnames, start, end, strand)]

			return(merged)
		})
	} else if (data_type == "tsr") {
		merged_samples <- map(sample_list, function(sample_group) {

			select_samples <- experiment@experiment$TSRs[sample_group]

			tsr_consensus <- select_samples %>%
				purrr::reduce(c) %>%
				GenomicRanges::reduce(ignore.strand = FALSE)

			merged <- map(select_samples, function(x) {
					overlap <- findOverlapPairs(query = tsr_consensus, subject = x)
					overlap <- as.data.table(overlap)
			})
			merged <- rbindlist(merged)[,
				.(first.seqnames, first.start, first.end,
				first.strand, second.X.score)
			]

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

			return(merged)
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
			rbindlist(id = "sample")
		
		tss_gr <- makeGRangesFromDataFrame(tss_set, keep.extra.columns = TRUE)

		# Make GRanges of TSRs.
		tsr_set <- extract_counts(experiment, "tsr", tsr_name) %>%
			rbindlist(id = "tsr_sample")
		
		tsr_gr <- tsr_set
		setnames(tsr_gr, old = c("FID", "FHASH"), new = c("TSR_FID", "TSR_FHASH"))
		tsr_gr <- tsr_gr[, .(seqnames, start, end, strand, TSR_FID, TSR_FHASH)]
		tsr_gr <- makeGRangesFromDataFrame(tsr_gr, keep.extra.columns = TRUE)

		# Annotate TSS to overlapping TSR.
		overlapping <- as.data.table(findOverlapPairs(tss_gr, tsr_gr))
		setnames(
			overlapping,
			old = c(
				"first.X.seqnames", "first.X.start", "first.X.end",
				"first.X.strand", "second.TSR_FID", "first.sample",
				"second.X.TSR_FHASH"
			),
			new = c("seqnames", "start", "end", "strand", "TSR_FID", "sample", "TSR_FHASH")
		)
		overlapping <- overlapping[, .(seqnames, start, end, strand, TSR_FID, TSR_FHASH, sample)]

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
		setkeyv(tsr_clean, c("TSR_FID", "TSR_FHASH"))
		setkeyv(tss_set, c("TSR_FID", "TSR_FHASH"))
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
