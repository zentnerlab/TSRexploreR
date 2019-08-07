
#' Annotate TSSs
#'
#' Use ChIPseeker package to annotate TSSs.
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom ChIPseeker annotatePeak
#' @importFrom purrr map
#' @importFrom GenomicFeatures makeTxDbFromGFF
#'
#' @param experiment tsrexplorer object with TSS Granges
#' @param upstream Bases upstream of TSS
#' @param downstream Bases downstream of TSS
#' @param gtf GTF genomic annotation file
#'
#' @return Annotated TSS slot in tsrchitect object
#'
#' @export
#' @rdname tss_annotation-function

tss_annotation <- function(experiment, upstream = 1000, downstream = 100, gtf = NA_character_) {
	## Load GTF.
	genome.annotation <- makeTxDbFromGFF(gtf, format = "gtf")
		
	## Annotate TSSs.
	tss.annotated <- map(
		experiment@experiment,
		~ annotatePeak(.,
			tssRegion = c(-upstream, downstream),
			TxDb = genome.annotation,
			sameStrand = TRUE
		) %>% as_tibble(.name_repair="unique")
	)
	
	experiment@annotated$TSSs <- tss.annotated
	return(experiment)
}
