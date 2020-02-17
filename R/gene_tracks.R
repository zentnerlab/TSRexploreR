
#' Gene Tracks by Gene
#'
#' @description
#' Generate gene tracks in GViz by gene name.
#'
#' @import Gviz
#'
#' @param experiment tsrexplorer object
#' @param genome_annotation Genome annotation GTF/GFF file, or TxDb object
#' @param feature_name Name of gene or transcript to plot
#' @param feature_type Either 'gene' or 'transcript'
#' @param tss_samples Names of samples to plot TSSs for, or 'none'
#' @param tsr_samples Names of samples to plot TSRs for, or 'none'
#' @param threshold TSSs and TSRs below threshold are excluded from plotting
#' @param upstream bases upstream to extend gene or promoter track
#' @param downstream bases downstream to extend gene or promoter track
#' @param promoter_only Instead of plotting the entire gene, plot the promoter region
#' @param use_cpm Use CPM normalized reads or not
#' @param track_colors Either a single color value for all tracks, or a vector of colors
#' @param axis_scale Relative size scale for axis text and title
#' @param ymax Maximum value on Y axis for all TSS tracks
#'
#' @rdname gene_tracks-function
#' @export

gene_tracks <- function(
	experiment, genome_annotation, feature_name, feature_type = "gene",
	tss_samples = "all", tsr_samples = "all", threshold = 1,
	upstream = 250, downstream = 250, promoter_only = FALSE,
	use_cpm = FALSE, track_colors = "black", axis_scale = 0.25,
	ymax = NA
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

	## Get ranges for requested genes.
	feature_ranges <- feature_ranges[feature_name, ]

	## Grab selected samples.
        use_tss <- !tss_samples %in% c("", " ", "none") & !is.na(tss_samples)
        use_tsr <- !tsr_samples %in% c("", " ", "none") & !is.na(tsr_samples)

	if (use_tss) {
		selected_TSSs <- extract_counts(experiment, "tss", tss_samples, use_cpm)
	}
	if (use_tsr) {
		selected_TSRs <- extract_counts(expermient, "tsr", tsr_samples, use_cpm)
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
	if (length(track_colors) > 1) {
		track_colors <- unlist(map(track_colors, ~ rep(., 2)))
	} else {
		track_colors <- rep(track_colors, length(split_TSSs))	
	}
	names(track_colors) <- names(split_TSSs)

	# Genome annotation track.
	genome_track <- GeneRegionTrack(
		anno, name = "", shape = "arrow",
		col = NA, fill = "black", showId = TRUE
	)

	# Data tracks.
	tss_tracks <- split_TSSs %>%
		imap(function(gr, sample_name) {
			if (is.na(ymax)) {
				data_track <- DataTrack(
					gr, name = sample_name, cex.title = axis_scale,
					cex.axis = axis_scale, col.histogram = track_colors[sample_name],
					fill.histogram = track_colors[sample_name]
				)
			} else {
                                data_track <- DataTrack(
                                        gr, name = sample_name, cex.title = axis_scale,
                                        cex.axis = axis_scale, col.histogram = track_colors[sample_name],
                                        fill.histogram = track_colors[sample_name], ylim = c(0, ymax)
                                )
			}
			return(data_track)
		})

	# Combine and plot tracks.
	plotTracks(
		c(genome_track, tss_tracks),
		chromosome = seqnames(feature_ranges),
		from = start(feature_ranges),
		to = end(feature_ranges),
		background.title = "white",
		col.title = "black",
		col.axis = "black",
		type = "histogram",
		baseline = 0,
		col.baseline = "black"
	)
}
