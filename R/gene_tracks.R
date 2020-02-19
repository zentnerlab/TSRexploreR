
#' Gene Tracks by Gene
#'
#' @description
#' Generate gene tracks in GViz by gene name.
#'
#' @import Gviz
#' @importFrom stringr str_count
#'
#' @param experiment tsrexplorer object
#' @param genome_annotation Genome annotation GTF/GFF file, or TxDb object
#' @param feature_name Name of gene or transcript to plot
#' @param feature_type Either 'gene' or 'transcript'
#' @param samples Names of samples to plot.
#' Append sample names with 'TSS:' or 'TSR:' for TSS and TSR tracks respectively.
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
	experiment, genome_annotation, feature_name, feature_type = "gene",
	samples = "all", threshold = 1, upstream = 250, downstream = 250,
	promoter_only = FALSE, use_cpm = FALSE, tss_colors = "black",
	tsr_colors = "black", axis_scale = 0.25, ymax = NA
) {
	
	## Prepare genome annotation.
	if (is(genome_annotation, "character")) {
		anno <- makeTxDbFromGFF(genome_annotation)
	} else if (is(genome_annotation, "TxDb")) {
		anno <- genome_annotation
	}
	options(ucscChromosomeNames = FALSE)

	## Grab appropriate ranges.
	if (promoter_only) {
		if (feature_type == "transcript") {
			feature_ranges <- promoters(transcripts(anno), upstream = upstream, downstream = downstream)
		} else {
			feature_ranges <- promoters(genes(anno), upstream = upstream, downstream = downstream)
		}
	} else {
		if (feature_type == "transcript") {
			feature_ranges <- transcripts(anno)
		} else {
			feature_ranges <- genes(anno)
		}
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
				tss_granges <- as.data.table(x)[
					score >= threshold,
					.(seqnames, start, end, strand, score)
				]
				tss_granges <- makeGRangesFromDataFrame(tss_granges, keep.extra.columns = TRUE)
				return(tss_granges)
			})
	}

	if (use_tsr) {
		selected_TSRs <- selected_TSRs %>%
			map(function(x) {
				tsr_granges <- as.data.table(x)[
					score >= threshold,
					.(seqnames, start, end, strand, score)
				]
				tsr_granges <- makeGRangesFromDataFrame(tsr_granges, keep.extra.columns = TRUE)
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
		if (length(track_colors) > 1) {
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
					cex.axis = axis_scale, col.histogram = track_colors[sample_name],
					fill.histogram = tss_colors[sample_name]
				)
			} else {
                               	data_track <- DataTrack(
                                       	gr, name = sample_name, cex.title = axis_scale,
                                       	cex.axis = axis_scale, col.histogram = tss_colors[sample_name],
                                       	fill.histogram = track_colors[sample_name], ylim = c(0, ymax)
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
	tracks <- map(samples, function(x) {
		track_name <- str_replace(x, "^(TSS:|TSR:)", "")
		if (str_detect(x, "^TSS:")) {
			select_track <- c(
				tss_tracks[[str_c(track_name, ".pos")]],
				tss_tracks[[str_c(track_name, ".neg")]]
			)
		} else if (str_detect(x, "^TSR:")) {
			select_track <- tsr_tracks[[track_name]]
		}
		return(select_track)
	})
	tracks <- c(genome_track, tracks)

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
