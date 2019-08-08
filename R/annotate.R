
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

tss_annotation <- function(experiment, annotation_file, feature_type = c("gene", "transcript"), upstream = 1000, downstream = 100) {
	## Load GTF.
	genome_annotation <- makeTxDbFromGFF(annotation_file, format = "gtf")
		
	## Annotate TSSs.
	tss_annotated <- map(
		experiment@experiment$TSSs,
		~ annotatePeak(.,
			tssRegion = c(-upstream, downstream),
			TxDb = genome_annotation,
			sameStrand = TRUE,
			level = feature_type,
		) %>% as_tibble(.name_repair="unique")
	)
	
	experiment@annotated$TSSs <- tss_annotated
	return(experiment)
}

#' Annotate TSRs
#'
#' Use ChIPseeker to annotate TSRs
#'
#' @include tsrexplorer.R
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom ChIPseeker annotatePeak
#' @importFrom purrr map
#' @importFrom GenomicFeatures makeTxDbFromGFF
#'
#' @param experiment tsrexplorer object with TSR Granges
#' @param upstream Bases upstream of TSR
#' @param downstream Bases downstream of TSR
#' @param gtf GTF genomic annotation file
#'
#' @return Annotated TSR slot in tsrchitect object
#'
#' @export
#' @rdname tsr_annotation-function

tsr_annotation <- function(experiment, annotation_file, feature_type = c("gene", "transcript"), upstream = 1000, downstream = 100) {
        ## Load GTF.
        genome_annotation <- makeTxDbFromGFF(annotation_file, format = "gtf")

        ## Annotate TSSs.
        tss_annotated <- map(
                experiment@experiment$TSRs,
                ~ annotatePeak(.,
                        tssRegion = c(-upstream, downstream),
                        TxDb = genome_annotation,
                        sameStrand = TRUE,
			level = feature_type
                ) %>% as_tibble(.name_repair="unique")
        )

        experiment@annotated$TSRs <- tss_annotated
        return(experiment)
}
