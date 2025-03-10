
#' detect_damage
#
#' Estimate the level of cytoplasmic RNA loss in a cell, a fundamental
#' proxy for damage in single cell RNA sequencing, through comparison
#' to cells of the input data where RNA loss, ranging from 0 to 100 %, has
#' been simulated. The true and simulated cells are merged and processed
#' before the quality control metrics,
#'
#' * log(non-zero features)
#' * log(total counts)
#' * mitochondrial proportion
#' * ribosomal proportion
#' * MALAT1 expression
#'
#' are computed and compared through principal component analysis (PCA) to
#' generate a distance matrix. The top related cells, or nearest neighbours
#' (NNs), defined by the matrix are retrieved for each cell. For true cells,
#' the proportion of NNs that are artificial (pANNs), i.e. simulated cells,
#' are found for each level of loss simulated.The level of loss where the pANNs
#' is highest is used to assign a predicted level of loss to each true cell.
#' Damage labels are then assigned to true cells if the predicted level of
#' loss is greater than or equal to 40 % or a user specified threshold.
#'
#'
#' @param count_matrix Matrix or dgCMatrix containing the counts from
#'  single cell RNA sequencing data.
#' @param organism String specifying the organism of origin of the input
#'  data where there are two standard options,
#'
#'  * "Hsap"
#'  * "Mmus"
#'
#'  If a user wishes to use a non-standard organism they must input a list
#'  containing strings for the patterns to match mitochondrial and ribosomal
#'  genes of the organism. If available, nuclear-encoded genes that are likely
#'  retained in the nucleus, such as in nuclear speckles, must also
#'  be specified. An example for humans is below,
#'
#'  * organism = c(mito_pattern = "^MT-",
#'                 ribo_pattern = "^(RPS|RPL)",
#'                 nuclear <- c("NEAT1","XIST", "MALAT1")
#'
#' * Default is "Hsap"
#' @param filter_threshold Numeric between 0 and 1 specifying above what
#'   proportion of estimated cytoplasmic RNA loss a cell should be regarded
#'   as damaged. We recommend not going lower than 0.4.
#'
#'  * Default 0.4
#' @param verbose Boolean specifying whether messages and function progress
#'  should be displayed in the console.
#'
#'  * Default is TRUE.
#' @return
#' @export
#'
#' @examples
detect_damage <- function(
    count_matrix,
    organism = "Hsap",
    filter_threshold = 0.4,
    verbose = TRUE
){

  # Phase One: Generate artificial damaged cells  ----

  # Phase Two: Combine artificial cells with true cells  ----

  # Compute quality control metrics ----

  # Phase Three: Run principal component analysis


  return(damage_labels)
}
