#' Merge Samples
#'
#' @description
#' Merge TSSs or TSRs by group or select samples.
#'
#' @inheritParams common_params
#' @param data_type Either 'tss' or 'tsr'.
#' @param merge_group Column in sample sheet to merge by.
#' @param merge_list Named list of samples to merge.
#'   List names will be the new TSS/TSR name, and the list contents should
#'   be a character vector of the TSSs or TSRs to merge.
#' @param merge_group The name of the column in the sample sheet
#'   that has the factor levels to merge samples by.
#' @param merge_list If merge_group is set to 'FALSE',
#'   specify what samples to merge in list format.
#'   List names will be the new TSS/TSR name, and the list contents should
#'   be a character vector of the TSSs or TSRs to merge.
#' @param max_distance Merge TSRs within this distance.
#'
#' @details
#' This function will merge overlapping TSSs or TSRs from different samples
#'   using either the sample sheet, or a named list.
#' If 'merge_group' is specified, the new merged TSS/TSR set will be the
#'   factor level in the column, and all TSS/TSR sets sharing that factor level
#'   will be merged.
#' If 'merge_list' is specified instead,
#'   The new TSS/TSR set will be the name of the list element,
#'   and the samples to merge will be a character vector as the list element.
#'
#' 'merge_distance' is provided for TSRs, and will merge TSRs within a certain
#'   distance from other TSRs.
#'
#' @return TSRexploreR object containing merged TSSs or TSRs.
#'
#' @examples
#' data(TSSs)
#' sample_sheet <- data.frame(
#'   sample_name=sprintf("S288C_D_%s", seq_len(3)),
#'   file_1=NA, file_2=NA,
#'   condition="Diamide"
#' )
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#'
#' tsre <- TSSs %>%
#'   tsr_explorer(sample_sheet=sample_sheet, genome_assembly=assembly) %>%
#'   format_counts(data_type="tss")
#'
#' # Merge TSSs by condition column.
#' merge_samples(tsre, data_type="tss", merge_group="condition")
#'
#' # Merge TSRs by condition column.
#' tsre <- tss_clustering(tsre, threshold=3)
#' merge_samples(tsre, data_type="tsr", merge_group="condition")
#'
#' @export

merge_samples <- function(
  experiment,
  data_type=c("tss", "tsr"),
  threshold=NULL,
  sample_sheet=NULL,
  merge_group=NULL,
  merge_list=NULL,
  max_distance=NULL,
  genome_assembly=NULL
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr")
  )
  assert_that(
    is.null(threshold) ||
    (is.numeric(threshold) && threshold >= 0)
  )
  assert_that(
    is.null(sample_sheet) ||
    (is.character(sample_sheet) | is.data.frame(sample_sheet))
  )
  assert_that(is.null(merge_group) || is.string(merge_group))
  assert_that(
    is.null(merge_list) ||
    (is.list(merge_list) && has_attr(merge_list, "names"))
  )
  if (is.null(merge_group) & is.null(merge_list)) {
    stop("Either 'merge_group' or 'merge_list' must be specified")
  }
  if (!is.null(merge_group) & !is.null(merge_list)) {
    stop("Either 'merge_group' or 'merge_list' must be specified")
  }

  ## If group is specified, prepare sample list.
  merge_type <- case_when(
    !is.null(merge_group) ~ "sample_sheet",
    !is.null(merge_list) ~ "sample_list"
  )

  if (merge_type == "sample_sheet") {
    keep <- c("sample_name", merge_group)
    merge_list <- experiment@meta_data$sample_sheet[, ..keep] %>%
      split(by=merge_group, keep.by=FALSE) %>%
      map(`[[`, "sample_name")
  }

  ## Merge feature sets.
  merged_samples <- map(merge_list, function(x) {

    # Bind ranges.
    samples <- experiment %>%
      extract_counts(data_type, x) %>%
      rbindlist(idcol="samples")

    # Merge ranges.
    norm_status <- "normalized_score" %in% colnames(samples)

    samples <- switch(
      data_type,
      "tss"=.merge_overlapping_TSSs(samples, norm_status),
      "tsr"=.merge_overlapping_TSRs(
        samples, norm_status, max_distance,
        experiment, genome_assembly
      )
    )

    # Convert to GRanges and ensure proper sorting.
    samples <- sort(as_granges(samples))

    return(samples)
    
  })

  ## Add merged GRanges to experiment slot.
  if (data_type == "tss") {
    experiment@experiment$TSSs <- c(experiment@experiment$TSSs, merged_samples)
  } else if (data_type == "tsr") {
    experiment@experiment$TSRs <- c(experiment@experiment$TSRs, merged_samples)
  }

  ## Add merged GRanges to count slot.
  merged_samples <- map(merged_samples, as.data.table)
  if (data_type == "tss") {
    experiment@counts$TSSs$raw <- c(experiment@counts$TSSs$raw, merged_samples)
  } else if (data_type == "tsr") {
    experiment@counts$TSRs$raw <- c(experiment@counts$TSRs$raw, merged_samples)
  }

  ## Add sample info to sample sheet.
  if (is.null(experiment@meta_data$sample_sheet)) {
    setDT(sample_sheet)
    experiment@meta_data$sample_sheet <- sample_sheet
  }

  return(experiment)
}

#' Merge Overlapping TSSs
#'
#' @param samples Sample data
#' @param norm_status Whether there is a normalized score column.

.merge_overlapping_TSSs <- function(
  samples,
  norm_status
) {

  if (norm_status) {
    samples <- samples[,
      .(score=sum(score), normalized_score=sum(normalied_score)),
      by=.(seqnames, start, end, strand)
    ]
  } else {
    samples <- samples[,
      .(score=sum(score)),
      by=.(seqnames, start, end, strand)
    ]
  }

  samples[, FHASH := str_c(seqnames, start, end, strand, sep=":")]

  return(samples)
}

#' Merge Overlapping TSRs
#'
#' @inheritParams common_params
#' @param samples Sample data.
#' @param norm_status Whether there is a normalized score column.
#' @param max_distance Maximum distance of TSRs to be merged.

.merge_overlapping_TSRs <- function(
  samples,
  norm_status,
  max_distance,
  experiment,
  genome_assembly
) {

  ## Open genome assembly.
  genome_assembly <- .prepare_assembly(genome_assembly, experiment)

  ## Create GRanges.
  samples_granges <- as_granges(samples)    

  ## Expand ranges if required and remove out of bounds.
  if (!is.null(max_distance)) {
    samples_granges <- stretch(samples_granges, max_distance * 2)

    # Add chromosome lengths to GRanges.
    assembly_type <- case_when(
      is(genome_assembly, "BSgenome") ~ "bsgenome",
      is(genome_assembly, "FaFile") ~ "fafile"
    )
  
    chrm_lengths <- switch(
      assembly_type,
      "fafile"=Rsamtools::seqinfo(genome_assembly),
      "bsgenome"=GenomeInfoDb::seqinfo(genome_assembly)
    )
  
    chrm_lengths <- chrm_lengths[seqlevels(samples_granges)]
    seqlengths(samples_granges) <- seqlengths(chrm_lengths)

    # Trim out-of-bounds ranges.
    samples_granges <- GenomicRanges::trim(samples_granges)
  }

  ## Create consensus ranges.
  cranges <- reduce_ranges_directed(samples_granges)
  cranges <- as.data.table(cranges, key=c("seqnames", "strand", "start", "end"))
  cranges[, group_ID := .I]

  ## Merge consensus ranges into original ranges.
  if (norm_status) {
    samples <- samples[, .(samples, seqnames, start, end, strand, score, normalized_score)]
    setkey(samples, seqnames, strand, start, end)
    overlaps <- foverlaps(cranges, samples)
  } else {
    samples <- samples[, .(samples, seqnames, start, end, strand, score)]
    setkey(samples, seqnames, strand, start, end)
    overlaps <- foverlaps(cranges, samples)
  }

  ## Aggregate scores.
  if (norm_status) {
    overlaps <- overlaps[, .(
      seqnames=unique(seqnames), start=min(start), end=max(end),
      strand=unique(strand), score=sum(score), normalized_score=sum(normalized_score)
    ), by=group_ID]
  } else {
    overlaps <- overlaps[, .(
      seqnames=unique(seqnames), start=min(start), end=max(end),
      strand=unique(strand), score=sum(score)
    ), by=group_ID]
  }
  overlaps[, group_ID := NULL]
  overlaps[, FHASH := str_c(seqnames, start, end, strand, sep=":")]

  return(overlaps)
}
