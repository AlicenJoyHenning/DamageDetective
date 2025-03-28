# tidyr pivot_longer variables
utils::globalVariables(c(
  "Features", "New_Features", "New_MitoProp", "New_RiboProp",
  "Original_Features", "Original_MitoProp", "Original_RiboProp", "Ribo. prop"
))


check_simulate_inputs <- function(
    count_matrix,
    damage_proportion,
    beta_shape_parameters,
    target_damage,
    distribution_steepness,
    damage_distribution,
    ribosome_penalty,
    organism
){
  # Check that count matrix is given
  if (is.null(count_matrix)) stop("Please provide 'count_matrix' input.")
  if (!inherits(count_matrix, "matrix") & !inherits(count_matrix, "CsparseMatrix")) {
    stop("Please ensure 'count_matrix' is a matrix or a sparse matrix (dgCMatrix).")
  }

  # Ensure user adjustments to default parameters are executable
  if (is.null(damage_proportion)) stop("Please provide 'damage_proportion'.")
  if (!is.numeric(damage_proportion) || damage_proportion < 0 || damage_proportion > 1) {
    stop("Please ensure 'damage_proportion' is a numeric between 0 and 1.")
  }
  if (!is.null(beta_shape_parameters) &
      length(beta_shape_parameters) != 2) {
    stop("Please ensure 'beta_shape_parameters' is a numeric vector of length 2.")
  }
  if (!is.numeric(target_damage) || length(target_damage) != 2 || target_damage[1] < 0 || target_damage[2] > 1 || target_damage[1] >= target_damage[2]) {
    stop("Please ensure 'target_damage' is a numeric vector of length 2, with values between 0 and 1, and the first value is less than the second.")
  }
  if (!distribution_steepness %in% c("shallow", "moderate", "steep")) {
    stop("Please ensure 'distribution_steepness' is one of 'shallow', 'moderate', or 'steep'.")
  }
  if (!damage_distribution %in% c("right_skewed", "left_skewed", "symmetric")) {
    stop("Please ensure 'damage_distribution' is one of 'right_skewed', 'left_skewed', or 'symmetric'.")
  }
  if (!is.numeric(ribosome_penalty) || ribosome_penalty < 0 || ribosome_penalty > 1) {
    stop("Please ensure 'ribosome_penalty' is a numeric between 0 and 1.")
  }
  if (!organism %in% c("Hsap", "Mmus") & length(organism) != 3) {
    stop("Please ensure 'organism' is one of 'Hsap' or 'Mmus', see documentation for non-standard organisms.")
  }
}

check_detect_inputs <- function(
    filter_threshold,
    mito_quantile,
    damage_levels,
    count_matrix,
    kN
){
  if (!is.numeric(filter_threshold) ||
      filter_threshold > 1 ||
      filter_threshold  < 0) {
    stop("Please ensure 'filter_threshold' is a numeric between 0 and 1.")
  }

  if (!is.numeric(mito_quantile) ||
      mito_quantile > 1 ||
      mito_quantile < 0) {
    stop("Please ensure 'mito_quantile' is a numeric between 0 and 1.")
  }

  if (!(damage_levels %in% c(3, 5, 7))) {
    if (!is.list(damage_levels) ||
        !all(grepl("^pANN_\\d+$", names(damage_levels))) ||
        !all(sapply(damage_levels, function(x) is.numeric(x) && length(x) == 2))) {
      stop("Please ensure `damage_levels` is of the correct format.")
    }
  }

  if (!is.null(kN)){
    if (!is.numeric(kN) ||
        kN > dim(count_matrix)[2]) {
      stop("Please ensure 'kN' is not greater than the number of cells.")
    }
  }
}

check_penalty_inputs <- function(
    mito_quantile,
    penalty_range,
    penalty_step,
    max_penalty_trials,
    stability_limit,
    return_output
){
  if (!is.numeric(mito_quantile) ||
      mito_quantile > 1 ||
      mito_quantile  < 0) {
    stop("Please ensure 'mito_quantile' is a numeric between 0 and 1.")
  }

  if (length(penalty_range) != 2) {
    stop("Please ensure 'penalty_range' is a numerical vector of length 2.")
  }

  if (penalty_range[[1]] > penalty_range[[2]] ||
      penalty_range[[1]] < 0 ||
      penalty_range[[2]] > 1){
    stop("Please ensure 'penalty_range' provides values between 0 and 1.")
  }

  if (!is.numeric(penalty_step) ||
      penalty_step < 0 ||
      penalty_step > 1){
    stop("Please ensure 'penalty_step' is a numeric between 0 and 1.")
  }

  if (!is.numeric(max_penalty_trials) ||
      max_penalty_trials > 1000){
    stop("Please ensure 'max_penalty_trials' lies within a reasonable range.")
  }

  if (!is.numeric(stability_limit) ||
      stability_limit > 1000){
    stop("Please ensure 'stability_limit' lies within a reasonable range.")
  }

  if (!return_output %in% c("penalty", "full")){
    stop("Please ensure 'return_output' is one of 'penalty' or 'full'.")
  }

}


#' Retrieve genes corresponding to the organism of interest
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
#' @return List containing the indices of the count matrix corresponding to
#'  mitochondrial, non-mitochondrial, and ribosomal gene sets.
#' @keywords internal
get_organism_indices <- function(
    count_matrix,
    organism
  ){
  if (organism == "Hsap") {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^(RPS|RPL)"
    nuclear <- c("FIRRE", "NEAT1","XIST", "MALAT1", "MIAT", "MEG3", "KCNQ1OT1", "HOXA11-AS", "FTX")
    MALAT1 <- "MALAT1"
  }

  if (organism == "Mmus") {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^(rps|rpl)"
    nuclear <- c("Firre", "Neat1","Xist", "Malat1", "Miat", "Meg3", "Kcnq1ot1", "Hoxa11-as", "Ftx")
    MALAT1 <- "Malat1"
  }

  # Allow for user specification for non-standard organism
  if (!organism %in% c("Hsap", "Mmus")) {
    mito_pattern <- organism$mito_pattern
    ribo_pattern <- organism$ribo_pattern
    nuclear <- organism$nuclear
    MALAT1 <- organism$MALAT1
  }

  # Isolate gene set indices (consistent across cells, not subsetting the matrix)
  mito_idx <- grep(mito_pattern, rownames(count_matrix), ignore.case = FALSE)
  nucl_idx <- which(rownames(count_matrix) %in% nuclear)
  mito_idx <- c(mito_idx, nucl_idx)
  non_mito_idx <- setdiff(seq_len(nrow(count_matrix)), mito_idx)
  ribo_idx <- grep(ribo_pattern, rownames(count_matrix), ignore.case = FALSE)

  return(list(
    mito_idx = mito_idx,
    ribo_idx = ribo_idx,
    non_mito_idx = non_mito_idx,
    mito_pattern =  mito_pattern,
    ribo_pattern = ribo_pattern,
    MALAT1 = MALAT1
  ))

}
