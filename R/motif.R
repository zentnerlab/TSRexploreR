#' Retrieve Sequences Near TSSs
#'
#' @description
#' Retrieve sequences surrounding TSSs for further plotting
#'
#' @include TSRexplore.R
#'
#' @importFrom tools file_ext
#'
#' @inheritParams common_params
#' @param distance Bases to add on each side of eacg TSS
#'
#' @details
#' This function will retrieve the genomic sequence surrounding TSSs for later use in
#'   plotting sequence color maps or sequence logos.
#'
#' 'genome_assembly' must be a valid genome assembly in either fasta or BSgenome format.
#' fasta formatted genome assemblies should have the file extension '.fasta' or '.fa'.
#' BSgenome assemblies are precompiled Bioconductor libraries for common organisms.
#'
#' 'distance' controls the length upstream and downstream of the TSS
#'   from which the sequence will be retrieved.
#'
#' A set of functions to control data structure for plotting are included.
#' 'threshold' will define the minimum number of reads a TSS or TSR
#'  must have to be considered.
#' 'dominant' specifies whether only the dominant TSS or TSR is considered 
#'   from the 'mark_dominant' function.
#' For TSSs this can be either dominant per TSR or gene, and for TSRs
#'   it is just the dominant TSR per gene.
#' 'data_conditions' allows for the advanced filtering, ordering, and grouping
#'   of data.
#'
#' @return DataFrame of sequences surrounding TSSs.
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#' seqs <- tss_sequences(tsre_exp, genome_assembly=assembly)
#'
#' @seealso
#' \code{\link{plot_sequence_logo}} to make sequence logos.
#' \code{\link{plot_sequence_colormap}} to make sequence color maps.

.tss_sequences <- function(
  experiment,
  samples="all",
  genome_assembly=NULL,
  threshold=NULL,
  use_normalized=FALSE,
  distance=10,
  dominant=FALSE,
  data_conditions=NULL
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(
    is.null(genome_assembly) || is.character(genome_assembly) ||
    is(genome_assembly, "BSgenome")
  )
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.count(distance))
  assert_that(is.flag(dominant))
  assert_that(is.null(data_conditions) || is.list(data_conditions))

  ## Open genome assembly.
  genome_assembly <- .prepare_assembly(genome_assembly, experiment)

  ## Get selected samples.
  select_samples <- experiment %>%
    extract_counts("tss", samples, use_normalized) %>%
    preliminary_filter(dominant, threshold)

  ## Apply data conditions.
  select_samples <- condition_data(select_samples, data_conditions)

  ## Prepare table for sequence retrieval.
  select_samples <- rbindlist(select_samples, idcol="sample")
  select_samples[, tss := start]
  select_samples <- as_granges(select_samples)

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

  chrm_lengths <- chrm_lengths[seqlevels(select_samples)]
  seqlengths(select_samples) <- seqlengths(chrm_lengths)

  ## Expand GRanges and remove out of bound.
  select_samples <- stretch(select_samples, distance * 2)

  out_of_bounds <- GenomicRanges:::get_out_of_bound_index(select_samples)
  if (length(out_of_bounds) > 0) {
    select_samples <- select_samples[-out_of_bounds]
  }

  ## Retrieve sequences.
  seqs <- switch(
    assembly_type,
    "bsgenome"=BSgenome::getSeq(genome_assembly, select_samples),
    "fafile"=Rsamtools::getSeq(genome_assembly, select_samples)
  )

  seqs <- seqs %>%
    as.data.table %>%
    {cbind(as.data.table(select_samples), .)}
      
  setnames(seqs, old="x", new="sequence")

  ## Order samples if required.
  if (!all(samples == "all")) {
    seqs[, sample := factor(sample, levels=samples)]
  }
  
  return(seqs)
}

#' Generate Sequence Logo
#'
#' @description
#' Create a sequence logo for the sequences around TSSs.
#'
#' @import ggseqlogo
#' @importFrom Biostrings consensusMatrix
#' @importFrom cowplot plot_grid
#'
#' @inheritParams common_params
#' @param distance Bases to add on each side of eacg TSS
#' @param font_size Font size for plots
#' @param base_colors Colors for each base.
#' @param ... Arguments passed to ggseqlogo.
#'
#' @details
#' This plotting function uses the ggseqlogo library to make sequence logos
#'   from the sequences retrieved by the 'tss_sequences' function.
#' Sequence logos show the enrichment of bases with certain positional biases
#'   in a centered set of sequences.
#' This is particularly important for TSS analysis since literature has shown
#'   strong base preferences spanning TSSs and surrounding sequences.
#'
#' @return ggplot2 object with sequence logo
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#' seqs <- tss_sequences(tsre_exp, genome_assembly=assembly)
#' plot_sequence_logo(seqs)
#'
#' @seealso
#' \code{\link{tss_sequences}} to get the surrounding sequences.
#' \code{\link{plot_sequence_colormap}} for a sequence color map plot.
#'
#' @rdname plot_sequence_logo-function
#' @export

plot_sequence_logo <- function(
  experiment,
  samples="all",
  genome_assembly=NULL,
  threshold=NULL,
  use_normalized=FALSE,
  distance=10,
  dominant=FALSE,
  data_conditions=NULL,
  ncol=1,
  font_size=6,
  base_colors=c(
    A="#109649", C="#255C99",
    G="#F7B32C", T="#D62839"
  ),
  ...
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(
    is.null(genome_assembly) || is.character(genome_assembly) ||
    is(genome_assembly, "BSgenome")
  )
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.count(distance))
  assert_that(is.flag(dominant))
  assert_that(is.null(data_conditions) || is.list(data_conditions))
  assert_that(is.count(ncol))
  assert_that(is.numeric(font_size) && font_size > 0)
  assert_that(
    is.character(base_colors) &&
    all(c("A", "T", "G", "C") %in% names(base_colors))
  )

  ## Get sequences.
  tss_sequences <- .tss_sequences(
    experiment,
    samples,
    genome_assembly,
    threshold,
    use_normalized,
    distance,
    dominant,
    data_conditions
  )

  ## Store status of data conditions.
  grouping_status <- case_when(
    !is.null(data_conditions$quantiling) ~ "row_quantile",
    !is.null(data_conditions$grouping) ~ "row_groups",
    TRUE ~ "none"
  )

  ## Grab sequences from input.
  if (grouping_status == "none") {
    sequences <- tss_sequences %>%
      split(.$sample) %>%
      map(function(x) {x[["sequence"]]})
  } else {
    setnames(tss_sequences, old=grouping_status, new="grouping")
    sequences <- tss_sequences %>%
      as.data.table %>%
      split(.$grouping) %>%
      map(function(x) {
        split(x, x$sample) %>%
        map(function(y) {y[["sequence"]]})
      })
  }

  ## Create viridis color scheme for bases.
  viridis_bases <- make_col_scheme(
    chars=c("A", "C", "G", "T"),
    groups=c("A", "C", "G", "T"),
    cols=base_colors[match(
      c("A", "C", "G", "T"),
      names(base_colors)
    )]
  )

  ## Make sequence logo.
  if (grouping_status == "none") {
    p <- ggseqlogo(sequences, ncol=ncol, ...) +
      theme(text=element_text(size=font_size)) #+
      #scale_x_continuous(
      # breaks=c(1, distance, distance + 1, (distance * 2) + 1),
      # labels=c(-distance, -1, +1, distance + 1)
      #)
  } else {
    p <- rev(sequences) %>%
      map(function(x) {
        ggseqlogo(x, ncol=ncol, ...) +
          theme(text=element_text(size=font_size)) #+
          #scale_x_continuous(
          # breaks=c(1, distance, distance + 1, (distance * 2) + 1),
          # labels=c(-distance, -1, +1, distance + 1)
          #)
      })

    p <- plot_grid(plotlist=p, labels=rev(names(sequences)), ncol=1)
  }

  return(p)
}

#' Plot Sequence Colormap
#'
#' Make a sequence colormap for the sequences around TSSs.
#'
#' @importFrom dplyr bind_cols
#' @importFrom stringr str_length
#'
#' @inheritParams common_params
#' @param base_colors Named vector specifying colors for each base
#' @param distance Bases to add on each side of eacg TSS
#' @param font_size Size of text for plots
#' @param ... Arguments passed to geom_tile
#'
#' @details
#' This plotting function generates a ggplot2 base color map surrounding TSSs.
#' Base color maps represent each base surrounding a TSS as a different color.
#' Since the base composition for every TSS can be seen in one plot, it's a good
#'   companion figure to sequence logos.
#'
#' The color of each base is set using the 'base_colors' argument.
#' The argument input should be a named vector, with the base as the name,
#' and the desired color of the base as the vector element.
#'
#' @return ggplot2 object of sequence colormap
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type="tss")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
#' seqs <- tss_sequences(tsre_exp, genome_assembly=assembly)
#' plot_sequence_colormap(seqs)
#'
#' @seealso
#' \code{\link{tss_sequences}} to get the surrounding sequence.
#' \code{\link{plot_sequence_logo}} to plot a sequence logo.
#'
#' @rdname plot_sequence_colormap-function
#' @export

plot_sequence_colormap <- function(
  experiment,
  samples="all",
  genome_assembly=NULL,
  threshold=NULL,
  use_normalized=FALSE,
  distance=10,
  dominant=FALSE,
  data_conditions=NULL,
  ncol=1,
  base_colors=c(
    "A"="#109649", "C"="#255C99",
    "G"="#F7B32C", "T"="#D62839"
  ),
  font_size=6,
  rasterize=FALSE,
  raster_dpi=150,
  ...
) {

  ## Check inputs.
  assert_that(is(experiment, "tsr_explorer"))
  assert_that(is.character(samples))
  assert_that(
    is.null(genome_assembly) || is.character(genome_assembly) ||
    is(genome_assembly, "BSgenome")
  )
  assert_that(is.null(threshold) || (is.numeric(threshold) && threshold >= 0))
  assert_that(is.count(distance))
  assert_that(is.flag(dominant))
  assert_that(is.null(data_conditions) || is.list(data_conditions))
  assert_that(is.count(ncol))
  assert_that(is.numeric(font_size) && font_size > 0)
  assert_that(
    is.character(base_colors) &&
    all(c("A", "T", "G", "C") %in% names(base_colors))
  )
  assert_that(is.flag(rasterize))
  assert_that(is.count(raster_dpi))

  ## Get sequences.
  tss_sequences <- .tss_sequences(
    experiment,
    samples,
    genome_assembly,
    threshold,
    use_normalized,
    distance,
    dominant,
    data_conditions
  )

  ## Store status of data conditions.
  grouping_status <- case_when(
    !is.null(data_conditions$quantiling) ~ "row_quantile",
    !is.null(data_conditions$grouping) ~ "row_groups",
    TRUE ~ "none"
  )

  ## Split sequences into columns
  split_seqs <- tss_sequences[, tstrsplit(sequence, split="")]

  setnames(
    split_seqs,
    old=sprintf("V%s", seq(1, (distance * 2) + 1)),
    new=as.character(c(seq(-distance, -1), seq(1, distance + 1)))
  )
  seq_data <- cbind(tss_sequences, split_seqs)

  ## Get order of TSSs for plotting.
  if (!is.null(data_conditions$ordering))  {
    seq_data[, FHASH := fct_reorder(FHASH, row_order)]
  }

  ## Format data for plotting.
  long_data <- melt(
      seq_data,
      measure.vars=as.character(c(seq(-distance, -1), seq(1, distance + 1))),
      variable.name="position", value.name="base"
    )

  long_data[,
    c("position", "base") := list(
      position=as.numeric(position),
      base=factor(base, levels=c("A", "C", "G", "T"))
    )
  ]
    
  ## Plot sequence colormap
  p <- ggplot(long_data, aes(x=.data$position, y=.data$FHASH))

  if (rasterize) {
    p <- p + rasterize(
      geom_tile(aes(fill=.data$base, color=.data$base), ...),
      dpi=raster_dpi
    )
  } else {
    p <- p + geom_tile(aes(fill=.data$base, color=.data$base), ...)
  }

  p <- p +
    scale_fill_manual(values=base_colors) +
    scale_color_manual(values=base_colors) +
    theme_minimal() +
    theme(
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      legend.title=element_blank(),
      axis.title.x=element_text(margin=margin(t=15)),
      panel.grid=element_blank(),
      text=element_text(size=font_size)
    ) +
    scale_x_continuous(
      breaks=c(1, distance, distance + 1, (distance * 2) + 1),
      labels=c(-distance, -1, 1, distance + 1)
    )

  if (grouping_status == "none") {
    p <- p + facet_wrap(. ~ sample, scales="free", ncol=ncol)
  } else {
    setnames(long_data, old=grouping_status, new="grouping")
    p <- p + facet_wrap(grouping ~ sample, scales="free", ncol=ncol)
  }

  return(p)
}
