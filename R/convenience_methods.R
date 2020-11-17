#' Add TSSs
#'
#' Convenience function to add TSSs to a TSRexploreR object.
#'
#' @include TSRexplore.R
#'
#' @param TSRexploreR_object TSRexploreR object
#'
#' @rdname tss_experiment-generic
#'
#' @export

setGeneric("tss_experiment", function(TSRexploreR_object) standardGeneric("tss_experiment"))
setGeneric("tss_experiment<-", function(TSRexploreR_object, value) standardGeneric("tss_experiment<-"))

#' @rdname tss_experiment-generic

setMethod("tss_experiment", signature(TSRexploreR_object="tsr_explorer"),
  function(TSRexploreR_object) {
    TSRexploreR_object@experiment$TSSs
  }
)

setMethod("tss_experiment<-", signature(TSRexploreR_object="tsr_explorer"),
  function(TSRexploreR_object, value) {
    TSRexploreR_object@experiment$TSSs <- value
  }
)

#' Add TSRs
#'
#' Convenience function to add TSRs to a TSRexploreR object.
#'
#' @include TSRexplore.R
#'
#' @param TSRexploreR_object TSRexploreR object
#'
#' @rdname tsr_experiment-generic
#'
#' @export

setGeneric("tsr_experiment", function(TSRexploreR_object) standardGeneric("tsr_experiment"))
setGeneric("tsr_experiment<-", function(TSRexploreR_object, value) standardGeneric("tsr_experiment<-"))

#' @rdname tsr_experiment-generic

setMethod("tsr_experiment", signature(TSRexploreR_object="tsr_explorer"),
  function(TSRexploreR_object) {
    TSRexploreR_object@experiment$TSRs
  }
)

setMethod("tsr_experiment<-", signature(TSRexploreR_object="tsr_explorer"),
  function(TSRexploreR_object, value) {
    TSRexploreR_object@experiment$TSRs <- value
  }
)

#' Extract Counts
#'
#' @description
#' Extract counts from a TSRexploreR object.
#'
#' @param experiment TSRexploreR object
#' @param data_type Whether to extract TSS or TSR counts
#' @param samples Names of samples from which to extract counts
#' @param use_normalized Whether CPM values should be returned
#'
#' @rdname extract_counts-function
#' @export

extract_counts <- function(experiment, data_type, samples, use_normalized=FALSE) {

  ## Extract appropriate TSS or TSR samples.
  if (data_type == "tss") {
    if (all(samples == "all")) samples <- names(experiment@counts$TSSs$raw)
    selected_samples <- experiment@counts$TSSs$raw[samples]
  } else if (data_type == "tsr") {
    if (all(samples == "all")) samples <- names(experiment@counts$TSRs$raw)
    selected_samples <- experiment@counts$TSRs$raw[samples]
  } else if (data_type == "tss_features") {
    if (all(samples == "all")) samples <- names(experiment@counts$TSS_features$raw)
    selected_samples <- experiment@counts$TSS_features$raw[samples]
  } else if (data_type == "tsr_features") {
    if (all(samples == "all")) samples <- names(experiment@counts$TSR_features$raw)
    selected_samples <- experiment@counts$TSR_features$raw[samples]
  }

  ## Return counts as copies so you don't unintentionally overwrite the TSRexploreR object values.
  return_samples <- map(selected_samples, copy)

  ## Return CPM score if requested.
  if (use_normalized) {
    walk(return_samples, function(x) {
      x[, score := normalized_score]
      x[, normalized_score := NULL]
      return(x)
    })
  }

  ## Return the ranges and counts as output of function.
  return(return_samples)
  
}

#' Extract Count Matrices
#'
#' Extract normalized count matrices.
#'
#' @param experiment TSRexploreR object
#' @param data_type Whether to extract 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param samples Samples to extract
#'
#' @rdname extract_matrix-function
#' @export

extract_matrix <- function(experiment, data_type, samples) {

  if (data_type == "tss") {
    selected_samples <- experiment@counts$TSSs$matrix
  } else if (data_type == "tsr") {
    selected_samples <- experiment@counts$TSRs$matrix
  } else if (data_type == "tss_features") {
    selected_samples <- experiment@counts$TSS_features$matrix
  } else if (data_type == "tsr_features") {
    selected_samples <- experiment@counts$TSR_features$matrix
  }

  if (samples == "all") samples <- colnames(selected_samples)
  selected_samples <- selected_samples[, samples]

  return(selected_samples)
}

#' Extract Differential Gene Sets
#'
#' Extract differential expression analysis results.
#'
#' @param experiment TSRexploreR object
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param de_comparisons The comparison sets to extract
#'
#' @rdname extract_de-function
#' @export

extract_de <- function(experiment, data_type, de_comparisons="all") {
        
  if (data_type == "tss") {
    if (de_comparisons == "all") {
      de_samples <- experiment@diff_features$TSSs$results
    } else {
      de_samples <- experiment@diff_features$TSSs$results[de_comparisons]
    }
  } else if (data_type == "tsr") {
    if (de_comparisons == "all") {
      de_samples <- experiment@diff_features$TSRs$results
    } else {
      de_samples <- experiment@diff_features$TSRs$results[de_comparisons]
    }
  } else if (data_type == "tss_features") {
    if (de_comparisons == "all") {
      de_samples <- experiment@diff_features$TSS_features$results
    } else {
      de_samples <- experiment@diff_features$TSS_features$results[de_comparisons]
    }
  } else if (data_type == "tsr_features") {
    if (de_comparisons == "all") {
      de_samples <- experiment@diff_features$TSR_features$results
    } else {
      de_samples <- experiment@diff_features$TSR_features$results[de_comparisons]
    }
  }

  return(de_samples)
}

#' Set slot data.
#'
#' Set the content of a slot.
#'
#' @param tsre_obj TSRexploreR object
#' @param new_data Data to be added to the slot
#' @param fill_slot either 'counts' or 'diff_features'
#' @param data_type Either 'tss', 'tsr', 'tss_features', or 'tsr_features'
#' @param subslot Either 'raw', 'matrix', or 'results'
#'
#' @rdname set_slot-function
#' @export

set_count_slot <- function(
  tsre_obj,
  new_data,
  fill_slot,
  data_type,
  subslot
) {

  ## Check inputs
  if (!is(tsre_obj, "tsr_explorer")) stop("tsre_obj must be a tsr explorer object")
  data_type <- match.arg(
    str_to_lower(data_type),
    c("tss", "tsr", "tss_features", "tsr_features", "tss_diff", "tsr_diff")
  )
  subslot <- match.arg(str_to_lower(subslot), c("raw", "matrix", "results"))

  ## Convert data type to their slot name.
  data_type <- switch(
    data_type,
    "tss"="TSSs",
    "tsr"="TSRs",
    "tss_diff"="TSSs",
    "tsr_diff"="TSRs",
    "tss_features"="TSS_features",
    "tsr_features"="TSR_features"
  )

  ## Fill the slot
  slot(tsre_obj, fill_slot)[[data_type]][[subslot]] <- new_data
  return(tsre_obj)
}
  
