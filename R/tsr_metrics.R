
#' TSR Metrics
#'
#' Calculate TSR Metrics
#'
#' @param experiment tsrexplorer object
#'
#' @rdname tsr_metrics-function
#' @export

tsr_metrics <- function(experiment) {

	## Grab samples from tsrexplorer object.
	select_samples <- experiment %>%
		extract_counts("tss", "all") %>%
		bind_rows(.id = "sample")

	keys <- c("sample", "TSR_FID")
	setkeyv(select_samples, keys)

	## Calculate shape index.
	si <- shape_index(select_samples)
	setkeyv(si, keys)
	select_samples <- merge(select_samples, si, all.x = TRUE)

	## Calculate peak concentration.
	pc <- peak_concentration(select_samples)
	setkeyv(pc, keys)
	select_samples <- merge(select_samples, pc, all.x = TRUE)

	## Calculate peak balance.
	pb <- peak_balance(select_samples)
	setkeyv(pb, keys)
	select_samples <- merge(select_samples, pb, all.x = TRUE)

	## Add metrics back to tsrexplorer object.
	experiment@counts$TSSs$raw <- select_samples
	return(experiment)
}

#' Shape Index
#'
#' Calculate shape index
#'
#' @param tss_table data.table of TSSs
#'
#' @rdname shape_index-function
#' @export

shape_index <- function(tss_table) {
	
	## Calculate shape index.
        si_results <- tss_table[
                !is.na(TSR_FID) & score > 1,
                .(shape_index = ((score / tsr_score) * log2(score / tsr_score))),
                by = .(sample, TSR_FID)
        ][,
                .(shape_index = 2 + sum(shape_index)), by = .(sample, TSR_FID)
        ]

        ## Classify based on shape index.
        si_results[, shape_class := ifelse(shape_index < -1, "broad", "peaked")]
	
	## Return results.
	return(copy(si_results)) 
}

#' Peak Concentration
#'
#' Calculate peak concentration
#'
#' @param tss_table data.table of TSSs
#'
#' @rdname peak_concentration-function
#' @export

peak_concentration <- function(tss_table) {
	
	## Calculate peak concentration.
	pc_results <- tss_table[
		!is.na(TSR_FID) & score > 1,
		.(peak_concentration = log2((max(score) / sum(score)) * max(score))),
		by = .(sample, TSR_FID)
	]

	## Return peak concentration results.
	return(copy(pc_results))
}

#' Peak Balance
#'
#' Calculate peak balance
#'
#' @param tss_table data.table of TSSs
#'
#' @rdname peak_balance-function
#' @export

peak_balance <- function(tss_table) {

	## Get TSS Position relative to TSR midpoints.
	tss_position <- tss_table[
		!is.na(TSR_FID) & score > 1,
		.(score, tsr_score, tss_pos = ifelse(
			strand == "+",
			start - median(range(start)),
			median(range(start)) - start
		)),
		by = .(sample, TSR_FID)
	]

	## Calculate peak balance.
	pb_results <- tss_position[,
		.(peak_balance = sum((score / sum(score)) * tss_pos)),
		by = .(sample, TSR_FID)
	]

	## Return peak balance values.
	return(copy(pb_results))

}
