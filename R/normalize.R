#' TMM Normalization
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
	counts <- experiment@experiment %>%
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
		spread(key = sample, value = score, fill = 0) %>%
		# Convert tibble to data frame, and turn the TSS name into rownames.
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
