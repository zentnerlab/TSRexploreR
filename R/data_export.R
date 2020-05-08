
#' Export TSSs
#'
#' Export TSSs to table or bedgraph
#'
#' @param experiment tsrexplorer object
#' @param samples Samples to export
#' @param file_type either 'bedgraph' or 'table'
#' @param out_dir Output directory for files
#'
#' @rdname tss_export-function
#'
#' @export

tss_export <- function(
	experiment,
	samples = "all",
	file_type = "bedgraph",
	out_dir = NA
) {

	## Retrieve samples.
	export_samples <- experiment@counts$TSSs$raw[samples]

	if (file_type == "bedgraph") {
		iwalk(export_samples, function(x, y) {
			x <- sort(makeGRangesFromDataFrame(x, keep.extra.columns = TRUE))
			
			pos_data <- x[strand(x) == "+"]
			pos_file <- file.path(
				ifelse(is.na(out_dir), getwd(), out_dir),
				str_c(y, "_pos.bedgraph")
			)
			export(pos_data, pos_file, "bedgraph")

			min_data <- x[strand(x) == "-"]
			min_file <- file.path(
				ifelse(is.na(out_dir), getwd(), out_dir),
				str_c(y, "_min.bedgraph")
			)
			export(min_data, min_file, "bedgraph")
		})
	} else if (data_type == "table") {
		iwalk(export_samples, function(x, y) {
			export_file <- file.path(
				ifelse(is.na(out_dir), getwd(), out_dir),
				str_c(y, "_TSSs.tsv")
			)
			fwrite(
				x, export_file, sep = "\t", col.names = TRUE,
				row.names = FALSE, quote = FALSE
			)
		})
	}
}
