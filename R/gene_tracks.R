#' Gene Tracks
#'
#' @description
#' Generate gene tracks in GViz by gene name.
#'
#' @importFrom Gviz GeneRegionTrack DataTrack AnnotationTrack plotTracks
#'   GenomeAxisTrack
#' @importFrom stringr str_count
#' @importFrom GenomicFeatures genes transcripts promoters
#'
#' @inheritParams common_params
#' @param feature_name Name of gene or transcript to plot.
#' @param feature_type Either 'gene' or 'transcript'
#' @param samples Names of samples to plot. Append sample names with 'TSS:' or 'TSR:' 
#' for TSS and TSR tracks, respectively.
#' @param upstream Bases upstream to extend gene or promoter track.
#' @param downstream Bases downstream to extend gene or promoter track.
#' @param promoter_only Instead of plotting the entire gene, plot only the promoter region.
#' @param tss_colors Either a single color value for all TSS tracks, or a vector of colors.
#' @param tsr_colors Either a single color value for all TSR tracks, or a vector of colors.
#' @param axis_scale Relative size scale for axis text and title.
#' @param ymax Maximum value on y-axis for all TSS tracks.
#' @param anno_pos Genome annotation and axis track position.
#'   Either 'top' or 'bottom'.
#'
#' @rdname gene_tracks-function
#' @export

gene_tracks <- function(
  experiment,
  feature_name,
  genome_annotation=NULL,
  feature_type="gene",
  samples="all",
  threshold=NULL,
  upstream=250,
  downstream=250,
  promoter_only=FALSE,
  use_normalized=FALSE,
  tss_colors="black",
  tsr_colors="black",
  axis_scale=0.25,
  ymax=NA,
  anno_pos="top"
) {

  ## Input checks.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(
    is.null(genome_annotation) || is.character(genome_annotation) ||
    is(genome_annotation, "TxDb")
  )
  assert_that(is.string(feature_name))
  feature_type <- match.arg(str_to_lower(feature_type), c("gene", "transcript"))
  assert_that(is.character(samples))
  names(samples) <- match.arg(
    str_to_lower(names(samples)),
    c("tss", "tsr"), several.ok=TRUE
  )
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.count(upstream))
  assert_that(is.count(downstream))
  assert_that(is.flag(promoter_only))
  assert_that(is.flag(use_normalized))
  assert_that(is.character(tss_colors))
  assert_that(is.character(tsr_colors))
  assert_that(is.numeric(axis_scale) && axis_scale > 0)
  assert_that(is.na(ymax) || is.numeric(ymax))
  anno_pos <- match.arg(
    str_to_lower(anno_pos),
    c("top", "bottom")
  )

  ## Prepare genome annotation.
  anno <- .prepare_annotation(genome_annotation, experiment)
  
  options(ucscChromosomeNames=FALSE)

  ## Grab appropriate ranges.
  ftype <- case_when(
    promoter_only & feature_type == "transcript" ~ "transcript_promoter",
    promoter_only & feature_type == "gene" ~ "gene_promoter",
    !promoter_only & feature_type == "transcript" ~ "transcript",
    !promoter_only & feature_type == "gene" ~ "gene"
  )

  feature_ranges <- switch(ftype,
    "transcript_promoter"=promoters(transcripts(anno), upstream=upstream, downstream=downstream),
    "gene_promoter"=promoters(genes(anno), upstream=upstream, downstream=downstream),
    "transcript"=transcripts(anno),
    "gene"=genes(anno)
  )


  ## If not promoter only, expand ranges.
  if (!promoter_only) {
    feature_ranges <- feature_ranges %>%
      anchor_5p %>%
      stretch(downstream) %>%
      anchor_3p %>%
      stretch(upstream)
  }

  ## Get ranges for requested gene.
  feature_ranges <- switch(feature_type,
    "transcript"=feature_ranges[feature_ranges$tx_name == feature_name, ],
    "gene"=feature_ranges[feature_ranges$gene_id == feature_name, ]
  )

  ## Get selected samples and convert to granges.
  use_tss <- any(names(samples) == "tss")
  use_tsr <- any(names(samples) == "tsr")

  keep_cols <- c("seqnames", "start", "end", "strand", "score")

  if (use_tss) {
    tss_samples <- samples[names(samples) == "tss"]
    selected_TSSs <- experiment %>%
      extract_counts("tss", tss_samples, use_normalized) %>%
      preliminary_filter(FALSE, threshold)
    selected_TSSs <- map(selected_TSSs, function(x) {
      x <- x[, ..keep_cols]
      x <- as_granges(x)
      return(x)
    })
  }
  if (use_tsr) {
    tsr_samples <- samples[names(samples) == "tsr"]
    selected_TSRs <- experiment %>%
      extract_counts("tsr", tsr_samples, use_normalized) %>%
      preliminary_filter(FALSE, threshold)
    selected_TSRs <- map(selected_TSRs, function(x) {
      x <- x[, ..keep_cols]
      x <- as_granges(x)
      return(x)
    })

  }

  ## Split positive and negative strands for TSSs.
  if (use_tss) {
    split_TSSs <- selected_TSSs %>%
      map(function(x) {
        pos_ranges <- plyranges::filter(x, strand == "+")
        neg_ranges <- plyranges::filter(x, strand == "-")
        split_ranges <- list(pos=pos_ranges, neg=neg_ranges)
        return(split_ranges)
      })
    split_TSSs <- unlist(split_TSSs, recursive=FALSE)
  }

  ## Build gene tracks.

  # Set track colors.
  if (use_tss) {
    if (length(tss_colors) > 1) {
      tss_colors <- rep(tss_colors, each=2)
    } else {
      tss_colors <- rep(tss_colors, length(split_TSSs))
    }
    names(tss_colors) <- names(split_TSSs)
  }

  if (use_tsr) {
    if (length(tsr_colors) == 1) {
      tsr_colors <- rep(tsr_colors, length(selected_TSRs))
    }
    names(tsr_colors) <- names(selected_TSRs)
  }

  # Genome annotation track.
  genome_track <- GeneRegionTrack(
    anno, name="", shape="arrow", col=NA, fill="black",
    showId=TRUE, cex.group=axis_scale
  )

  # Data tracks.
  if (use_tss) {
    tss_tracks <- imap(split_TSSs, function(gr, sample_name) {
      track_args <- list(
        gr, name=sample_name, cex.title=axis_scale, cex.axis=axis_scale,
        col.histogram=tss_colors[sample_name],
        fill.histogram=tss_colors[sample_name]
      )
      if (!is.na(ymax)) track_args <- c(track_args, list(ylim=c(0, ymax)))

      data_track <- do.call(DataTrack, track_args)
      return(data_track)
    })
  }

  if (use_tsr) {
    tsr_tracks <- imap(selected_TSRs, function(gr, sample_name) {
      anno_track <- AnnotationTrack(
        gr, name=sample_name, fill=tsr_colors[sample_name],
        cex.title=axis_scale, col=NA 
      )
    })
  }

  # Combine and plot tracks.
  tracks <- imap(samples, function(x, y) {
    track <- list()
    if (y == "tss") {
      track[[str_c(x, ".pos")]] <- tss_tracks[[str_c(x, ".pos")]]
      track[[str_c(x, ".neg")]] <- tss_tracks[[str_c(x, ".neg")]]
    } else if (y == "tsr") {
      track[[x]] <- tsr_tracks[[x]]
    }
    return(track)
  })

  tracks <- purrr::flatten(tracks)

  if (anno_pos == "top") {
    tracks <- c(
      list("axis_track"=GenomeAxisTrack()),
      list("genome_track"=genome_track),
      tracks
    )
  } else if (anno_pos == "bottom") {
    tracks <- c(
      tracks,
      list("genome_track"=genome_track),
      list("axis_track"=GenomeAxisTrack())
    )
  }

  plotTracks(
    tracks,
    chromosome=seqnames(feature_ranges),
    from=start(feature_ranges),
    to=end(feature_ranges),
    background.title="white",
    col.title="black",
    col.axis="black",
    type="histogram",
    baseline=0,
    col.baseline="black",
    title.width=2
  )
}
