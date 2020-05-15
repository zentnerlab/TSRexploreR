
#' Average Plots
#'
#' @description
#' Generate average plots of TSSs or TSRs
#'
#' @import tibble
#' @import data.table
#' @import ggplot2
#' @importFrom magrittr %>% extract
#' @importFrom purrr pmap
#'
#' @param experiment tsr_explorer object with annotated TSSs
#' @param samples Either 'all' to plot all samples, or a vector of sample names
#' @param data_type Make average plot for TSS or TSR data
#' @param consider_score Should the TSS or TSR score be considered, or just the unique location.
#' @param upstream Bases upstream of plot center
#' @param downstream Bases downstream of plot center
#' @param threshold threshold value for TSSs
#' @param ncol Number of columns to use when plotting data when quantiles not set
#' @param use_cpm Whether to use the CPM normalized or raw counts
#' @param dominant Consider only dominant TSS or TSR
#' @param data_conditions Data conditioning filters
#' @param ... Arguments passed to geom_density
#'
#' @details
#' This plotting function generates an average plot of TSS or TSR signal
#'   relative to annotated TSSs.
#' The plot is returned as a ggplot2 object.
#'
#' By default only the TSS or TSR position is considered, effectively giving every
#'   TSS or TSR a score of 1.
#' If 'consider_score' is set to TRUE, the score of the TSS or TSR will be considered when
#'   making the average plot, which gives more weight to stronger TSSs or TSRs.
#'
#' The region around the annotated TSS used for plotting is controlled by
#'   'upstream' and 'downstream', which should be positive integers.
#'
#' A set of functions to control data structure for plotting are included.
#' 'use_cpm' will use the CPM normalized scores, which only matters if
#'   'consider_score' is TRUE.
#' 'threshold' will define the minimum number of reads a TSS or TSR
#'  must have to be considered.
#' 'dominant' specifies whether only the dominant TSS or TSR is considered 
#'   from the 'mark_dominant' function.
#' For TSSs this can be either dominant per TSR or gene, and for TSRs
#'   it is just the dominant TSR per gene.
#' 'data_conditions' allows for the advanced filtering, ordering, and grouping
#'   of data.
#'   
#'
#' @return ggplot2 object of average plot
#'
#' @examples
#' TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
#' TSSs <- readRDS(TSSs)
#' tsre_exp <- tsr_explorer(TSSs)
#' tsre_exp <- format_counts(tsre_exp, data_type = "tss")
#' annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
#' tsre_exp <- annotate_features(
#'   tsre_exp, annotation_data = annotation,
#'   data_type = "tss", feature_type = "transcript"
#' )
#' plot_average(tsre_exp, data_type = "tss")
#'
#' @seealso
#' \code{\link{annotate_features}} to annotate the TSSs or TSRs.
#'
#' @rdname plot_average-function
#' @export

plot_average <- function(
        experiment,
        data_type = c("tss", "tsr"),
        samples = "all",
        consider_score = FALSE,
        upstream = 1000,
        downstream = 1000,
        threshold = 1,
        ncol = 1,
        use_cpm = FALSE,
        dominant = FALSE,
        data_conditions = NA,
        ...
) {

	## Check inputs.
	if (!is(experiment, "tsr_explorer")) stop("experiment must be a tsrexplorer object")

        if (!is(data_type, "character")) stop("data_type must be a character")
        if (length(data_type) > 1) stop("data_type must be a character")
        data_type <- str_to_lower(data_type)
        if (!data_type %in% c("tss", "tsr")) {
                stop("data_type must be either 'tss' or 'tsr'")
        }

	if (!is(samples, "character")) stop("samples must be a character")

	if (!is(consider_score, "logical") | !is(use_cpm, "logical") | !is(dominant, "logical")) {
		stop("consider_score, use_cpm, and/or dominant must be logical")
	}

        if (!is(upstream, "numeric") | !is(downstream, "numeric")) {
                stop("upstream and downstream must be positive integers")
        }
        if (upstream %% 1 != 0 | downstream %% 1 != 0) {
                stop("upstream and downstream must be positive integers")
        }
        if (upstream < 0 | downstream < 0) stop("upstream and downstream must be positive integers")

	if (
                !is.na(threshold) && !is(threshold, "numeric") ||
                threshold %% 1 != 0 || threshold < 1
        ) {
                stop("threshold must be a positive integer greater than or equal to 1")
        }

	if (!is(ncol, "numeric") || ncol %% 1 != 0 || ncol < 1) {
		stop("ncol must be a positive integer")
	}

	if (is.na(data.conditions) && !is(data_conditions, "list")) {
		stop("data_conditions should be a list of values")
	}

        ## Assign color type.
        if (data_type == "tss") {
                color_type <- "#431352"
        } else if (data_type == "tsr") {
                color_type <- "#34698c"
        }

        ## Pull data out of appropriate slot.
        sample_data <- extract_counts(experiment, data_type, samples, use_cpm)

        ## Preliminary data preparation.
        if (dominant | !is.na(threshold)) {
                sample_data <- preliminary_filter(sample_data, dominant, threshold)
        }

        sample_data <- map(sample_data, function(x)  {
                x <- x[dplyr::between(distanceToTSS, -upstream, downstream)]
                return(x)
        })

        ## Condition data.
        if (!is.na(data_conditions)) {
                sample_data <- do.call(group_data, c(list(signal_data = sample_data), data_conditions))
        }

        ## Update data if score is also considered instead of just unique position.
        sample_data <- rbindlist(sample_data, idcol = "sample")
        if (consider_score) sample_data <- sample_data[rep(seq_len(.N), score)]

	## Set sample order if required.
	if (!all(samples == "all")) {
		sample_data[, samples := factor(samples, levels = samples)]
	}

        ## Plot averages.
        groupings <- any(names(data_conditions) %in% c("quantile_by", "grouping"))

        p <- ggplot(sample_data, aes(distanceToTSS)) +
                geom_density(fill = color_type, color = color_type, ...) +
                labs(
                        x = "Position Relative to Annotated TSS",
                        y = "Density"
                ) +
                theme_bw()

        if (groupings) {
                p <- p + facet_grid(fct_rev(factor(grouping)) ~ sample)
        } else {
                p <- p + facet_wrap(~ sample, ncol = ncol)
        }
        return(p)
}
