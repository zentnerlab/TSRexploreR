
#' Group Data
#'
#' Group data for plotting
#'
#' @param signal_data TSS or TSR data
#' @param quantile_by Continuous metric to calculate quantiles
#' @param n_quantiles Number of quantiles to calculate for continuous metric
#' @param order_by Metric to order data by
#' @param order_direction Whether the values should be ordered in
#' 'ascending' or 'descending' order
#' @param filter_by Continuous variable to filter the data by
#' @param upper_filter Upper value to filter continuous data by
#' @param lower_filter Lower value to filter continuous data by
#'
#' @rdname group_data-function
#' @export

group_data <- function(
	signal_data,
	filter_by = NA, upper_filter = NA, lower_filter = NA,
	order_by = NA, order_direction = "descending",
	quantile_by = NA, n_quantiles = NA
) {
	
	## First filter the data if requested.	
	if (!is.na(filter_by)) {
		signal_data <- map(signal_data, function(x) {
			x <- x[!is.na(x[[filter_by]])]

			if (!is.na(upper_filter)) {
				x <- x[x[[filter_by]] <= upper_filter]
			}
			if (!is.na(lower_filter)) {
				x <- x[x[[filter_by]] >= lower_filter]
			}

			return(x)
		})
	}

	## Next order the data for plotting if requested.
	if (!is.na(order_by)) {
		signal_data <- map(signal_data, function(x) {
			x <- x[!is.na(x[[order_by]])]

			if (order_direction == "descending") {
				x[, plot_order := dense_rank(desc(x[[order_by]]))]
			} else if (order_direction == "ascending") {
				x[, plot_order := dense_rank(x[[order_by]])]
			}

			return(x)
		})
	}

	## Quantile the metric if requested.
	if (!is.na(quantile_by)) {
		signal_data <- map(signal_data, function(x) {
			x <- x[!is.na(x[[quantile_by]])]

			x[, ntile := ntile(x[[quantile_by]], n_quantiles)]

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
