#' Group Data
#'
#' Group data for plotting
#'
#' @param signal_data TSS or TSR data
#' @param quantile_by Continuous metric to calculate quantiles
#' @param n_quantiles Number of quantiles to calculate for continuous metric
#' @param quantile_samples Samples to use when setting quantiles
#' @param quantile_group Group these features and quantile based on aggregate mean of feature
#' @param order_by Metric to order data by
#' @param order_direction Whether the values should be ordered in
#' 'ascending' or 'descending' order
#' @param order_on Names of samples to order by
#' @param filters Logical string to subset/filter data by
#' @param grouping If quantiles not set split data by categorical variable
#'
#' @rdname group_data-function
#' @export

group_data <- function(
        signal_data, filters = NA,
        order_by = NA, order_direction = "descending",
        order_samples = NA, order_group = NA,
        quantile_by = NA, n_quantiles = NA,
	quantile_samples = NA, quantile_group = NA,
	quantile_direction = "descending",
        grouping = NA
) {
	
        ## First filter the data if requested.  
        if (!is.na(filters)) {
                signal_data <- map(signal_data, function(x) {
                        x <- subset(x, eval(parse(text = filters)))
                        return(x)
                })
        }

        ## Next order the data for plotting if requested.
        if (!is.na(order_by)) {
                order_dataset <- rbindlist(signal_data, id = "sample")[,
                        c("sample", ..order_by, ..order_group)
                ]

                order_dataset <- dcast(
                        order_dataset, as.formula(str_c("sample ~ ", order_group)),
                        value.var = order_by, fun.aggregate = mean,
                        fill = ifelse(order_direction == "descending", -Inf, Inf)
                )

                order_dataset <- melt(
                        order_dataset, id.vars = "sample",
                        variable.name = order_group, value.name = order_by
                )

                if (!is.na(order_samples)) {
                        order_dataset <- order_dataset[sample %in% order_samples]
                }

                order_dataset <- order_dataset[, .(order_by = mean(get(order_by))), by = order_group]
                setnames(order_dataset, old = "order_by", new = order_by)

                if (order_direction == "descending") {
                        order_dataset[, plot_order := dense_rank(desc(order_dataset[[order_by]]))]
                } else if (order_direction == "ascending") {
                        order_dataset[, plot_order := dense_rank(order_dataset[[order_by]])]
                }

                order_dataset[, c(order_by) := NULL]

                signal_data <- map(signal_data, function(x) {
                        x <- merge(x, order_dataset, by = order_group, all.x = TRUE)
                        return(x)
                })
        }

	## Quantile the metric if requested.
	if (!is.na(quantile_by)) {
		
		quantile_dataset <- rbindlist(signal_data, id = "sample")[,
			c("sample", ..quantile_by, ..quantile_group)
		]

		quantile_dataset <- dcast(
			quantile_dataset, as.formula(str_c("sample ~ ", quantile_group)),
			value.var = quantile_by, fun.aggregate = mean,
			fill = ifelse(quantile_direction == "descending", -Inf, Inf)
		)

		quantile_dataset <- melt(
			quantile_dataset, id.vars = "sample",
			variable.name = quantile_group, value.name = quantile_by
		)

		if (!is.na(quantile_samples)) {
			quantile_dataset <- quantile_dataset[sample %in% quantile_samples]
		}

		quantile_dataset <- quantile_dataset[, .(quantile_by = mean(get(quantile_by))), by = quantile_group]
		setnames(quantile_dataset, old = "quantile_by", new = quantile_by)
		
		if (quantile_direction == "descending") {
			quantile_dataset[, grouping := ntile(desc(quantile_dataset[[quantile_by]]), n_quantiles)]
		} else if (quantile_direction == "ascending") {
			quantile_dataset[, grouping := ntile(quantile_dataset[[quantile_by]], n_quantiles)]
		}

		quantile_dataset[, c(quantile_by) := NULL]

		signal_data <- map(signal_data, function(x) {
			x <- merge(x, quantile_dataset, by = quantile_group, all.x = TRUE)
			return(x)
		})
	}

	## If not using quantiles split by a categorical value.
	if (!is.na(grouping) & is.na(quantile_by)) {
		signal_data <- map(signal_data, function(x) {
			x <- x[!is.na(x[[grouping]])]

			x[, grouping := x[[grouping]]]

			return(x)
		})
	}

	## Return the signal data.
	return(signal_data)
}

#' Preliminary Filter
#'
#' Preliminary filter of data
#'
#' @param signal_data TSS or TSR data
#' @param dominant Whether to retain dominant TSS/TSR
#' @param threshold Wheher to apply a threshold to the scores
#'
#' @rdname preliminary_filter-function
#' @export

preliminary_filter <- function(signal_data, dominant, threshold) {
	
	## Retain only dominant features if requested.
	if (dominant) {
		signal_data <- map(signal_data, function(x) {
			x <- x[(dominant)]
			return(x)
		})
	}

	## Apply a threshold to score if requested.
	if (!is.na(threshold)) {
		signal_data <- map(signal_data, function(x) {
			x <- x[score >= threshold]
			return(x)
		})
	}

	## Return singal data.
	return(signal_data)
}
