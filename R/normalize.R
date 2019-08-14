#' TMM Normalize TSSs
#'
#' Using edgeR to TMM normalize TSSs
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom dplyr bind_rows mutate select
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom tidyr spread
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object
#'
#' @return tibble of TMM normalized read counts
#'
#' @export
#' @rdname tss_normalization-function

tss_normalization <- function(experiment) {
	raw_counts <- experiment@experiment$TSSs %>%
		# Make the TSS name a concatenation of the chromosome, start, end, and strand.
		map(
			~as_tibble(., .name_repair = "unique") %>%
			mutate(position = paste(seqnames, start, end, strand, sep="_")) %>%
			select(position, score)
		) %>%
		# Select the TSS name column and the score.
		# Turn the list of tibbles into one tibble with a column specify what tibble the row came from.
		bind_rows(.id = "sample") %>%
		# Turn samples into column and TSS names into rows.
		spread(key = sample, value = score, fill = 0)

	experiment@raw_counts$TSSs <- raw_counts

	counts <- raw_counts %>%
		# Convert to data frame and set position as row names.
		as.data.frame %>%
		column_to_rownames("position") %>%
		# Convert data frame to count matrix.
		as.matrix

	tmm_tbl <- counts %>%
		# Create edger object.
		DGEList %>%
		# TMM normalize read counts.
		calcNormFactors %>%
		# Extract TMM normalized read counts.
		cpm %>%
		# Turn TMM normalized read counts into tibble.
		as_tibble(rownames="TSS_position", .name_repair="unique")

	experiment@normalized_counts$TSSs <- tmm_tbl
	return(experiment)
}

#' TMM Normalize TSRs
#'
#' Using edgeR to TMM normalize TSRs
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom GenomicRanges GRangesList reduce findOverlaps makeGRangesFromDataFrame
#' @importFrom dplyr bind_rows mutate select group_by summarize
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom tidyr complete
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @param experiment tsrexplorer object with TSR GRanges
#'
#' @return tibble of TMM normalized read counts
#'
#' @export
#' @rdname tsr_normalization-function

tsr_normalization <- function(experiment) {
	## Merge overlapping TSRs to get consensus
	tsr_consensus <- experiment@experiment$TSRs %>%
		as("GRangesList") %>%
		unlist %>%
		reduce(ignore.strand=FALSE) %>%
		as_tibble(.name_repair = "unique") %>%
		mutate(names = paste(seqnames, start, end, strand, sep="_")) %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)
	
	names(tsr_consensus) <- sprintf("TSR_%.6d", 1:length(tsr_consensus))

	## Get overlapping TSRs with consensus
	overlapping <- map(
		names(experiment@experiment$TSRs),
		~findOverlaps(
			query = tsr_consensus,
			subject = makeGRangesFromDataFrame(experiment@experiment$TSRs[[.x]])
		)  %>%
			as_tibble(.name_repair = "unique") %>%
			mutate(
				nTAGs = experiment@experiment$TSRs[[.x]][subjectHits]$nTAGs,
				TSR_name = tsr_consensus[queryHits]$names
			) %>%
			select(-queryHits, -subjectHits) %>%
			group_by(TSR_name) %>%
			summarize(nTAGs = sum(nTAGs)) %>%
			complete(TSR_name = names(tsr_consensus), fill = list(nTAGs = 0))
			
	) %>% setNames(names(experiment@experiment$TSRs))

	## Create count matrix
	raw_counts <- overlapping %>%
		bind_rows(.id = "sample") %>%
		spread(key = sample, value = nTAGs)

	experiment@raw_counts$TSRs <- raw_counts

	count_matrix <- raw_counts %>%
		as.data.frame %>%
		column_to_rownames("TSR_name") %>%
		as.matrix

	## TMM normalize TSR matrix
	tmm_tbl <- count_matrix %>%
		DGEList %>%
		calcNormFactors %>%
		cpm %>%
		as_tibble(rownames="consensus_TSR", .name_repair="unique")

	experiment@normalized_counts$TSRs <- tmm_tbl
	return(experiment)
}
