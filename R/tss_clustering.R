
#' TSS Clustering
#'
#' Basic function to cluster TSSs
#'
#' @param experiment tsrexplorer object
#' @param threshold Minimum number of reads for a TSS to be considered
#' @param samples The samples to call TSRs for
#' @param max_distance The maximum distance to cluster TSSs
#'
#' @rdname tss_clustering-function
#' @export

tss_clustering <- function(
	experiment, samples = "all", threshold = 1, max_distance = 25
) {

	## Retrive samples.
	select_samples <- extract_counts(experiment, "tss", samples)
	
	select_samples <- select_samples %>%
		map(function(x) {
			x <- x[
				score >= threshold,
				.(seqnames, start, end, strand, score)
			]
			return(x)
		})

	## Convert samples to GRanges.
	select_samples <- map(
		select_samples,
		~ makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
	)

	## Call TSRs.
	clustered_TSSs <- select_samples %>%
		map(function(x) {

			# Cluster TSSs within 'max_distance'.
			clustered <- GenomicRanges::reduce(
				x, with.revmap = TRUE,
				min.gapwidth = max_distance + 1
			)

			# Get aggregate sum of scores.
			cluster_info <- aggregate(
				x, mcols(clustered)$revmap,
				score = sum(score),
				n_unique = lengths(score)
			)

			clustered$score <- cluster_info$score
			clustered$n_unique <- cluster_info$n_unique
			clustered$revmap <- NULL

			return(clustered)
		})

	## Add TSRs back to tsrexplorer object.
	clustered_TSSs <- map(clustered_TSSs, function(x) {
		x <- as.data.table(x)
		x[, FID := seq_len(nrow(x))]
		x[,
			FHASH := digest(str_c(seqnames, start, end, strand, collapse = "")),
			by = seq_len(nrow(x))
		]
		return(x)
	})

	TSR_granges <- map(clustered_TSSs, function(x) {
		TSR_granges <- x[, .(seqnames, start, end, strand, score)]
		TSR_granges <- makeGRangesFromDataFrame(TSR_granges, keep.extra.columns = TRUE)
		return(TSR_granges)
	})

	experiment@experiment$TSRs <- TSR_granges
	experiment@counts$TSRs$raw <- clustered_TSSs
	return(experiment)

}
