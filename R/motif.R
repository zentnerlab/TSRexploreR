
#' Retrieve Sequences Near TSSs
#'
#' @description
#' Retrieve sequences surrounding TSSs for further plotting
#'
#' @include tsrexplorer.R
#'
#' @importFrom tools file_ext
#' @importFrom Rsamtools FaFile
#'
#' @param experiment tsrexplorer object with TSS GRanges
#' @param samples Either "all" or a vector of names of samples to analyze
#' @param genome_assembly Genome assembly in fasta format or bioconductor BSgenome
#' @param threshold Keep only TSSs with at least this number of raw counts
#' @param distance Bases to add on each side of eacg TSS
#' @param dominant Whether to only consider dominant TSSs
#' @param data_conditions Condition the data (filter, quantile, and grouping available)
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
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "tsrexplorer")
#' seqs <- tss_sequences(tsre_exp, genome_assembly = assembly)
#'
#' @seealso
#' \code{\link{plot_sequence_logo}} to make sequence logos.
#' \code{\link{plot_sequence_colormap}} to make sequence color maps.
#'
#' @rdname tss_sequences-function
#' @export

tss_sequences <- function(
  experiment,
  samples = "all",
  genome_assembly,
  threshold = 1,
  distance = 10,
  dominant = FALSE,
  data_conditions = list(order_by = "score")
) {

  ## Check inputs.
  if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsrexplorer object")
  if (!is(samples, "character")) stop("samples must be a character")
  if (!is(genome_assembly, "character") & !is(genome_assembly, "BSgenome")) {
    stop("genome assembly must be a fasta file or BSgenome package")
  }
  if (is(genome_assembly, "character")) {
    extension <- file_ext(genome_assembly)
    if (!extension %in% c("fasta", "fa")) {
      stop("genome_assembly file extension must be '.fa' or '.fasta'")
    }
  }
  if (!is(threshold, "numeric") | !is(distance, "numeric")) {
    stop("threshold and/or distance must be positive integers")
  }
  if (threshold %% 1 != 0 | distance %% 1 != 0) {
    stop("threshold and/or distance must be positive integers")
  }
  if (threshold < 1 | distance < 1) {
    stop("threshold and/or distance must be positive integers")
  }
  if (!is(dominant, "logical")) stop("dominant must be logical")
  if (all(!is.na(data_conditions)) && !is(data_conditions, "list")) {
    stop("data_conditions must be a list of values")
  }

  ## Open genome assembly.
  assembly_type <- case_when(
    is(genome_assembly, "character") ~ "character",
    is(genome_assembly, "BSgenome") ~ "bsgenome"
  )

  genome_assembly <- switch(
    assembly_type,
    "character"=FaFile(genome_assembly),
    "bsgenome"=genome_assembly
  )

  ## Get selected samples.
  select_samples <- extract_counts(experiment, "tss", samples)

  ## Preliminary filtering of data.
  select_samples <- preliminary_filter(select_samples, dominant, threshold)

  ## Condition the data.
  if (all(!is.na(data_conditions))) {
    select_samples <- do.call(group_data, c(list(signal_data = select_samples), data_conditions))
  }

  ## Prepare table for sequence retrieval.
  select_samples <- rbindlist(select_samples, idcol = "sample")
  select_samples[, tss := start]

  ## Get sequences.
  seqs <- select_samples %>%
    as_granges %>%
    stretch(distance * 2) %>%
    {getSeq(genome_assembly, .)} %>%
    as.data.table %>%
    {cbind(select_samples, .)}
      
  setnames(seqs, old = "x", new = "sequence")

  ## Order samples if required.
  if (!all(samples == "all")) {
    seqs[, sample := factor(seqs, levels = samples)]
  }
  
  ## Generate and return DataFrame.
  groupings <- any(names(data_conditions) %in% c("quantile_by", "grouping"))
  
  keep_cols <- c(
    "sample", "FHASH", "grouping",
    "plot_order", "tss", "sequence", "score"
  )
  keep_cols <- keep_cols[keep_cols %in% colnames(seqs)]
  seqs <- seqs[, ..keep_cols]
  
  seqs <- DataFrame(seqs)
  metadata(seqs)$groupings <- groupings
  metadata(seqs)$threshold <- threshold
  metadata(seqs)$distance <- distance
  metadata(seqs)$dominant <- dominant

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
#' @param tss_sequences Sequences surrounding TSSs generated with tss_sequences
#' @param ncol Number of columns to plot if quantiles is not set
#' @param font_size Font size for plots
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
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "tsrexplorer")
#' seqs <- tss_sequences(tsre_exp, genome_assembly = assembly)
#' plot_sequence_logo(seqs)
#'
#' @seealso
#' \code{\link{tss_sequences}} to get the surrounding sequences.
#' \code{\link{plot_sequence_colormap}} for a sequence color map plot.
#'
#' @rdname plot_sequence_logo-function
#' @export

plot_sequence_logo <- function(
  tss_sequences,
  ncol = 1,
  font_size = 6
) {

  ## Check inputs.
  if (!is(tss_sequences, "DataFrame")) stop("tss_sequences must be a DataFrame")
  if (!is(ncol, "numeric") || ncol %% 1 != 0 || ncol < 1) {
    stop("ncol must be a positive integer")
  }
  if (!is(font_size, "numeric") || !font_size > 0) {
    stop("font_size must be a positive number")
  }

  ## Get some info used to pull sequencs.
  distance <- metadata(tss_sequences)$distance
  groupings <- metadata(tss_sequences)$groupings

  ## Grab sequences from input.
  if (!groupings) {
    sequences <- tss_sequences %>%
      as.data.table %>%
      split(.$sample) %>%
      map(function(x) {x[["sequence"]]})
  } else {
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
    chars = c("A", "C", "G", "T"),
    groups = c("A", "C", "G", "T"),
    cols = c("#431352", "#34698c", "#44b57b", "#fde540")
  )

  ## Make sequence logo.
  if (!groupings) {
    p <- ggseqlogo(sequences, ncol = ncol) +
      theme(text = element_text(size = font_size)) #+
      #scale_x_continuous(
      # breaks = c(1, distance, distance + 1, (distance * 2) + 1),
      # labels = c(-distance, -1, +1, distance + 1)
      #)
  } else {
    p <- rev(sequences) %>%
      map(function(x) {
        ggseqlogo(x, ncol = ncol) +
          theme(text = element_text(size = font_size)) #+
          #scale_x_continuous(
          # breaks = c(1, distance, distance + 1, (distance * 2) + 1),
          # labels = c(-distance, -1, +1, distance + 1)
          #)
      })

    p <- plot_grid(plotlist = p, labels = rev(names(sequences)), ncol = 1)
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
#' @param tss_sequences Sequences surrounding TSS generated with tss_sequences
#' @param ncol Number of columns to plot data if quantiles not specified
#' @param base_colors Named vector specifying colors for each base
#' @param text_size Size of text for plots
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
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' assembly <- system.file("extdata", "S288C_Assembly.fasta", package = "tsrexplorer")
#' seqs <- tss_sequences(tsre_exp, genome_assembly = assembly)
#' plot_sequence_colormap(seqs)
#'
#' @seealso
#' \code{\link{tss_sequences}} to get the surrounding sequence.
#' \code{\link{plot_sequence_logo}} to plot a sequence logo.
#'
#' @rdname plot_sequence_colormap-function
#' @export

plot_sequence_colormap <- function(
  tss_sequences,
  ncol = 1,
  base_colors = c(
    "A" = "#109649",
    "C" = "#255C99",
    "G" = "#F7B32C",
    "T" = "#D62839"
  ),
  text_size = 6,
  ...
) {
  ## Check inputs.
  if (!is(tss_sequences, "DataFrame")) stop("tss_sequences must be a DataFrame")
  if (!is(ncol, "numeric") || ncol %% 1 != 0 || ncol < 1) {
    stop("ncol must be a positive integer")
  }
  if (!is(text_size, "numeric") || !text_size > 0) {
    stop("font_size must be a positive number")
  }
  if (
    !is(base_colors, "character") || length(base_colors) < 4 ||
    is.null(names(base_colors)) ||
    !all(str_to_lower(names(base_colors)) %in% c("a", "t", "g", "c"))
  ) {
    stop("base_colors must be a named character vector")
  }

  ## Grab some information out of DataFrame.
  distance <- metadata(tss_sequences)$distance
  groupings <- metadata(tss_sequences)$groupings

  ## Start preparing data for plotting.
  seq_data <- as.data.table(tss_sequences)
  seq_data[, width := str_length(sequence)]

  ## Split sequences into columns
  split_seqs <- seq_data[, tstrsplit(sequence, split="")]

  #split_seqs <- seq_data[, as.data.table(str_split(sequence, "", simplify = TRUE))]
  setnames(
    split_seqs,
    old = sprintf("V%s", seq(1, (distance * 2) + 1)),
    new = as.character(c(seq(-distance, -1), seq(1, distance + 1)))
  )
  seq_data <- cbind(seq_data, split_seqs)

  ## Get order of TSSs for plotting.
  seq_data[, FHASH := fct_reorder(factor(FHASH), plot_order)]

  ## Format data for plotting.
  long_data <- seq_data %>%
    melt(
      measure.vars = as.character(c(seq(-distance, -1), seq(1, distance + 1))),
      variable.name = "position", value.name = "base"
    )

  long_data[,
    c("position", "base") := list(
      position = as.numeric(position),
      base = factor(base, levels = c("A", "C", "G", "T"))
    )
  ]
    
  ## Plot sequence colormap
  p <- ggplot(long_data, aes(x = .data$position, y = .data$FHASH)) +
    geom_tile(aes(fill = .data$base, color=.data$base)) +
    scale_fill_manual(values = base_colors) +
    scale_color_manual(values=base_colors) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      legend.title = element_blank(),
      axis.title.x = element_text(margin = margin(t = 15)),
      panel.grid = element_blank(),
      text = element_text(size = text_size)
    ) +
    scale_x_continuous(
      breaks = c(1, distance, distance + 1, (distance * 2) + 1),
      labels = c(-distance, -1, 1, distance + 1)
    )

  if (!groupings) {
    p <- p + facet_wrap(. ~ sample, scales = "free", ncol = ncol)
  } else {
    p <- p + facet_wrap(grouping ~ sample, scales = "free", ncol = ncol)
  }

  return(p)
}
