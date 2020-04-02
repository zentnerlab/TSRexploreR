
#' Sample Sheet
#'
#' Add sample sheet to tsrexplorer object
#'
#' @param experiment tsrexplorer object
#' @param sample_sheet Sample sheet as tab delimited file or data.frame
#' with columns: 'replicate_id', 'tss_name', 'tsr_name', and optionally 'group_id'
#'
#' @rdname add_sample_sheet-function
#' @export

add_sample_sheet <- function(experiment, sample_sheet) {

	## Import sample sheet if file.
	if (is(sample_sheet, "character")) {
		sample_sheet <- fread(sample_sheet, sep = "\t", header = TRUE)
	}

	## Convert sample sheet to data.table
	sample_sheet <- as.data.table(sample_sheet)
	
	## Add sample sheet to tsrexplorer object.
	experiment@meta_data$sample_sheet <- sample_sheet

	return(experiment)
}

#' Feature Groups
#'
#' Add gene/transcript group info to tsrexplorer object
#'
#' @param experiment tsrexplorer object
#' @param group_sheet Group sheet as tab delimited file or data.frame
#' with column 'feature', and other named columns containing feature groups
#' @param feature_type either 'gene' or 'transcript'
#'
#' @rdname add_feature_group-function
#' @export

add_feature_group <- function(experiment, group_sheet, feature_type) {
	
	## Import group sheet if file.
	if (is(group_sheet, "character")) {
		group_sheet <- read.delim(
			group_sheet, sep = "\t", header = TRUE,
			stringsAsFactors = FALSE
		)
	}

	## Convert group sheet to data.table
	group_sheet <- as.data.table(group_sheet)

	## Add group sheet to tsrexplorer object.
	experiment@meta_data$group_sheet <- group_sheet
	
	return(experiment)
}
