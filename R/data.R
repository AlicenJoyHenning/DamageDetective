#' Test Counts
#'
#' A sparse matrix of PBMC single-cell gene expression data for function
#' testing.
#'
#' This dataset consists of an abridged sparse matrix originally derived from
#' PBMC single-cell RNA sequencing data.
#' The data was processed and stored in a sparse matrix format.
#'
#' @format A sparse matrix with:
#' \itemize{
#'   \item 32 738 rows representing genes.
#'   \item 500 columns representing single cells.
#'   \item Stored expression values representing nonzero gene expression levels.
#' }
#'
#' @usage data(test_counts)
#' @name test_counts
#' @docType data
#' @source \doi{10.1038/s41591-020-0769-8}
"test_counts"

#' Penalty Plot
#'
#' A precomputed `ggplot` object visualizing the effects of different
#' ribosomal penalties on simulated damaged cell populations.
#'
#' This object consists of a combined patchwork plot with three panels:
#' one for each ribosomal penalty level (no penalty, medium penalty, heavy penalty),
#' generated using `DamageDetective::simulate_counts()`. These plots were created
#' to illustrate how ribosomal content penalization influences the separation
#' between healthy and damaged cells.
#'
#' @format A `gg` object created with the `ggplot2` and `patchwork` packages.
#'
#' @usage data(penalty_plot)
#' @name penalty_plot
#' @docType data
#' @seealso \code{\link[DamageDetective]{simulate_counts}}
"penalty_plot"

#' Selected Penalty
#'
#' Precomputed output of the penalty selection process where the full output
#' is returned in addition to the penalty selected.
#'
#' This object consists of a list with two items:
#' - A data frame containing the full penalty results including the global means
#' of each penalty tested.
#' - The value of the selected penalty
#'
#' @format A list containing a data frame and numeric.
#'
#' @usage data(selected_penalty)
#' @name selected_penalty
#' @docType data
#' @seealso \code{\link[DamageDetective]{select_penalty}}
"selected_penalty"

#' Filtered Output
#'
#' Abridged precomputed output of the detect damage process where filtering is
#' performed automatically.
#'
#' This object contains a numeric specifying the dimensions of the filtered
#' matrix generating.
#'
#' @format Numeric
#'
#' @usage data(filtered_output)
#' @name filtered_output
#' @docType data
#' @seealso \code{\link[DamageDetective]{detect_damage}}
"filtered_output"
