
#' Retrieve Sequences
#'
#' @inheritParams common_params
#' @param data_type Either 'tss', 'tsr', or 'shift'
#' @param fixed_size Whether all returned sequences should have the
#'   same fixed size starting from the center of the ranges
#' @param extend_upstream Distance to extend ranges upstream
#' @param extend_downstream Distance to extend ranges downstream
#' @param return_format Either 'table' or 'biostrings'
#' @param output_dir If set the results will be saved to this directory.
#'   If return format is set to 'biostrings' the sequences will be saved as
#'   FASTA files. If return format is a table it will be saved as a
#'   tab delimited table.
#'
#' @export

retrieve_seqs <- function(
  experiment,
  data_type=c("tss", "tsr", "shift"),
  samples="all",
  genome_assembly=NULL,
  threshold=NULL,
  dominant=FALSE,
  fixed_size=FALSE,
  extend_upstream=NULL,
  extend_downstream=NULL,
  return_format="biostrings",
  output_dir=NULL
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "shfit")
  )
  assert_that(is.character("samples"))
  assert_that(
    is.null(genome_assembly) ||
    is(genome_assembly, "BSgenome") ||
    is.character(genome_assembly)
  )
  assert_that(
    is.null(threshold) ||
    (is.numeric(threshold) && threshold >= 0)
  )
  assert_that(is.flag(dominant))
  assert_that(
    is.null(extend_upstream) ||
    is.count(extend_upstream)
  )
  assert_that(is.flag(fixed_size))
  assert_that(
    is.null(extend_downstream) ||
    is.count(extend_downstream)
  )
  return_format <- match.arg(
    str_to_lower(return_format),
    c("table", "biostrings")
  )
  assert_that(is.null(output_dir) || is.string(output_dir))

  ## Get selected samples.
  selected_samples <- experiment %>%
    extract_counts(data_type, samples) %>%
    preliminary_filter(dominant, threshold) %>%
    rbindlist(idcol="sample") %>%
    as_granges

  ## Prepare genome assembly.
  genome_assembly <- .prepare_assembly(genome_assembly, experiment)

  ## Add chromosome lengths to GRanges.
  assembly_type <- case_when(
    is(genome_assembly, "BSgenome") ~ "bsgenome",
    is(genome_assembly, "FaFile") ~ "fafile"
  )

  chrm_lengths <- switch(
    assembly_type,
    "fafile"=Rsamtools::seqinfo(genome_assembly),
    "bsgenome"=GenomeInfoDb::seqinfo(genome_assembly)
  )

  chrm_lengths <- chrm_lengths[seqlevels(selected_samples)]
  seqlengths(selected_samples) <- seqlengths(chrm_lengths)

  ## Reduce ranges to center point if requested.
  if (fixed_size) {
    selected_samples <- selected_samples %>%
      anchor_center %>%
      mutate(width=1)
  }

  ## Expand ranges downstream if requested.
  if (!is.null(extend_downstream)) {
    selected_samples <- selected_samples %>%
      anchor_5p %>%
      stretch(extend_downstream)
  }

  ## Expand ranges upstream if requested.
  if (!is.null(extend_upstream)) {
    selected_samples <- selected_samples %>%
      anchor_3p %>%
      stretch(extend_upstream)
  }

  ## Remove any ranges that are now out of bounds.
  out_of_bounds <- .out_of_bounds_index(selected_samples)
  if (length(out_of_bounds) > 0) {
    selected_samples <- selected_samples[-out_of_bounds]
  }

  ## Split the samples back into a list.
  selected_samples <- selected_samples %>%
    as.data.table %>%
    split(by="sample", keep.by=FALSE) %>%
    map(as_granges)

  ## Retrieve the sequences.
  seqs <- switch(
    assembly_type,
    "bsgenome"=map(selected_samples, ~BSgenome::getSeq(genome_assembly, .x)),
    "fafile"=map(selected_samples, ~Rsamtools::getSeq(genome_assembly, .x))
  )

  ## Add names back to biostrings.
  seqs <- map2(
    seqs, selected_samples,
    function(x, y) {
      names(x) <- y$FHASH
      return(x)
    }
  )

  ## Convert to table if required.
  if (return_format == "table") {
    seqs <- map2(
      seqs, selected_samples,
      function(x, y) {
        x <- as.data.frame(x)
        setDT(x, keep.rownames="FHASH")
        setnames(x, old="x", new="seq")

        y <- as.data.table(y)
        x <- merge(x, y, by="FHASH")

        return(x)
      }
    )
  }

  ## Save files if output directory is set.
  if (!is.null(output_dir)) {
    iwalk(seqs, function(x, y) {
      file_name <- switch(
        return_format,
        "biostrings"=file.path(output_dir, str_c(y, ".fasta")),
        "table"=file.path(output_dir, str_c(y, ".tsv"))
      )
      switch(
        return_format,
        "table"=fwrite(x, file_name, sep="\t"),
        "biostrings"=writeXStringSet(x, file_name, format="fasta")
      )
    })
  }

  ## Return the sequences if output directory is not set.
  if (is.null(output_dir)) {
    return(seqs)
  }

}

#' Get index of out of bounds ranges.
#'
#' @importFrom GenomeInfoDb isCircular
#'
#' @description
#' Taken from GenomicRanges:::get_out_of_bound_index
#'
#' @param x GenomicRanges object.

.out_of_bounds_index <- function(x) {
    if (length(x) == 0L)
        return(integer(0))
    x_seqnames_id <- as.integer(seqnames(x))
    x_seqlengths <- unname(seqlengths(x))
    seqlevel_is_circ <- unname(isCircular(x)) %in% TRUE
    seqlength_is_na <- is.na(x_seqlengths)
    seqlevel_has_bounds <- !(seqlevel_is_circ | seqlength_is_na)
    which(seqlevel_has_bounds[x_seqnames_id] &
          (start(x) < 1L | end(x) > x_seqlengths[x_seqnames_id]))
}
