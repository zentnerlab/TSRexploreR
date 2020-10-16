
#' Gene Tracks by Gene
#'
#' @description
#' Generate gene tracks in GViz by gene name.
#'
#' @importFrom Gviz GeneRegionTrack DataTrack AnnotationTrack plotTracks
#' @importFrom stringr str_count
#' @importFrom GenomicFeatures genes transcripts promoters
#'
#' @param experiment tsrexplorer object
#' @param genome_annotation Genome annotation GTF/GFF file, or TxDb object
#' @param feature_name Name of gene or transcript to plot
#' @param feature_type Either 'gene' or 'transcript'
#' @param samples Names of samples to plot.
#'   Append sample names with 'TSS:' or 'TSR:' for TSS and TSR tracks respectively.
#' @param threshold TSSs and TSRs below threshold are excluded from plotting
#' @param upstream bases upstream to extend gene or promoter track
#' @param downstream bases downstream to extend gene or promoter track
#' @param promoter_only Instead of plotting the entire gene, plot the promoter region
#' @param use_cpm Use CPM normalized reads or not
#' @param tss_colors Either a single color value for all TSS tracks, or a vector of colors
#' @param tsr_colors Either a single color value for all TSR tracks, or a vector of colors
#' @param axis_scale Relative size scale for axis text and title
#' @param ymax Maximum value on Y axis for all TSS tracks
#'
#' @rdname gene_tracks-function
#' @export

gene_tracks <- function(
  experiment,
  genome_annotation,
  feature_name,
  feature_type = "gene",
  samples = "all",
  threshold = 1,
  upstream = 250,
  downstream = 250,
  promoter_only = FALSE,
  use_cpm = FALSE,
  tss_colors = "black",
  tsr_colors = "black",
  axis_scale = 0.25,
  ymax = NA
) {

  ## Input checks.
  if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsr explorer object")
  if (!is(genome_annotation, "character") & !is(genome_annotation, "TxDb")) {
    stop("genome_annotation must be a file path or TxDb object")
  }
  if (!is(feature_name, "character") || length(feature_name) > 1) stop("feature_name must be a character")
  feature_type <- match.arg(str_to_lower(feature_type), c("gene", "transcript"))
  if (!is(samples, "character")) stop("samples must be a character")
  if (!all(str_detect(sample, "^TS[RS]:"))) stop("sample names must be prefixed with 'TSS:' or 'TSR:'")
  if (
    !is.na(threshold) && (!is(threshold, "numeric") ||
    threshold %% 1 != 0 || threshold < 1)
  ) {
    stop("threshold must be a positive integer")
  }
  if (!is(upstream, "numeric") | !is(downstream, "numeric")) {
    stop("upstream and downstream must be positive integers")
  }
  if (upstream %% 1 != 0 | downstream %% 1 != 0) {
    stop("upstream and downstream must be positive integers")
  }
  if (upstream < 0 | downstream < 0) stop("upstream and downstream must be positive integers")
  if (!is(promoter_only, "logical")) stop("promoter_only must be TRUE or FALSE")
  if (!is(use_cpm, "logical")) stop("use_cpm must be TRUE or FALSE")
  if (!is(tss_colors, "character")) stop("tss_colors must be a character vector")
  if (!is(tsr_colors, "character")) stop("tsr_colors must be a character vector")
  if (!is(axis_scale, "numeric") || length(axis_scale) > 1 || axis_scale <= 0) {
    stop("axis_scale must be a positive number")
  }
  if (!is.na(ymax) && (!is(ymax, "numeric") || length(ymax) > 1 || ymax <= 0)) {
    stop("ymax must be a positive number")
  }
  
  ## Prepare genome annotation.
  anno_type <- case_when(
    is(genome_annotation, "character") ~ "character",
    is(genome_annotation, "TxDb") ~ "txdb"
  )

  anno <- switch(
    anno_type,
    "character"=makeTxDbFromGFF(genome_annotation),
    "txdb"=genome_annotation
  )
  
  options(ucscChromosomeNames = FALSE)

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
      resize(., width = width(.) + downstream, "start") %>%
      resize(., width = width(.) + upstream,  "end")
  }

  ## Get ranges for requested gene.
  if (feature_type == "transcript") {
    feature_ranges <- feature_ranges[feature_ranges$tx_name == feature_name, ]
  } else if (feature_type == "gene") {
    feature_ranges <- feature_ranges[feature_ranges$gene_id == feature_name, ]
  }

  ## Grab selected samples.
  use_tss <- sum(str_count(samples, "^TSS:")) > 0
  use_tsr <- sum(str_count(samples, "^TSR:")) > 0

  if (use_tss) {
    tss_samples <- samples %>%
      keep(~ str_detect(., "^TSS:")) %>%
      str_replace("^TSS:", "")
    selected_TSSs <- extract_counts(experiment, "tss", tss_samples, use_cpm)
  }
  if (use_tsr) {
    tsr_samples <- samples %>%
      keep(~ str_detect(., "^TSR:")) %>%
      str_replace("^TSR:", "")
    selected_TSRs <- extract_counts(experiment, "tsr", tsr_samples, use_cpm)
  }

  ## Convert samples to GRanges.
  if (use_tss) {
    selected_TSSs <- selected_TSSs %>%
      map(function(x) {
        tss_granges <- x[
          score >= threshold,
          .(seqnames, start, end, strand, score)
        ]
        tss_granges <- as_granges(tss_granges)
        return(tss_granges)
      })
  }

  if (use_tsr) {
    selected_TSRs <- selected_TSRs %>%
      map(function(x) {
        tsr_granges <- x[
          score >= threshold,
          .(seqnames, start, end, strand, score)
        ]
        tsr_granges <- as_granges(tsr_granges)
        return(tsr_granges)
      })
  }

  ## Split positive and negative strands for TSSs.
  if (use_tss) {
    split_TSSs <- selected_TSSs %>%
      map(function(x) {
        pos_ranges <- x[strand(x) == "+",]
        neg_ranges <- x[strand(x) == "-",]
        split_ranges <- list(pos = pos_ranges, neg = neg_ranges)
        return(split_ranges)
      })
    split_TSSs <- unlist(split_TSSs, recursive = FALSE)
  }

  ## Build gene tracks.

  # Assign colors to tracks.
  if (use_tss) {
    if (length(tss_colors) > 1) {
      tss_colors <- unlist(map(tss_colors, ~ rep(., 2)))
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
    anno, name = "", shape = "arrow", col = NA, fill = "black",
    showId = TRUE, cex.group = axis_scale
  )

  # Data tracks.
  if (use_tss) {
    tss_tracks <- imap(split_TSSs, function(gr, sample_name) {
      if (is.na(ymax)) {
        data_track <- DataTrack(
          gr, name = sample_name, cex.title = axis_scale,
          cex.axis = axis_scale, col.histogram = tss_colors[sample_name],
          fill.histogram = tss_colors[sample_name]
        )
      } else {
        data_track <- DataTrack(
          gr, name = sample_name, cex.title = axis_scale,
          cex.axis = axis_scale, col.histogram = tss_colors[sample_name],
          fill.histogram = tss_colors[sample_name], ylim = c(0, ymax)
        )
      }
      return(data_track)
    })
  }

  if (use_tsr) {
    tsr_tracks <- imap(selected_TSRs, function(gr, sample_name) {
      anno_track <- AnnotationTrack(
        gr, name = sample_name, fill = tsr_colors[sample_name],
        cex.title = axis_scale, col = NA 
      )
    })
  }

  # Combine and plot tracks.
  tracks <- list("genome_track" = genome_track)
  for (samp in samples) {
    track_name <- str_replace(samp, "^(TSS:|TSR:)", "")
    if (str_detect(samp, "^TSS:")) {
      tracks[[str_c(track_name, ".pos")]] <- tss_tracks[[str_c(track_name, ".pos")]]
      tracks[[str_c(track_name, ".neg")]] <- tss_tracks[[str_c(track_name, ".neg")]]
    } else if (str_detect(samp, "^TSR:")) {
      tracks[[track_name]] <- tsr_tracks[[track_name]]
    }
  }

  plotTracks(
    tracks,
    chromosome = seqnames(feature_ranges),
    from = start(feature_ranges),
    to = end(feature_ranges),
    background.title = "white",
    col.title = "black",
    col.axis = "black",
    type = "histogram",
    baseline = 0,
    col.baseline = "black",
    title.width = 2
  )
}
