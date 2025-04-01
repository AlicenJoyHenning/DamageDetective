#' select_penalty
#'
#' Recommended prerequisite function to `detect_damage()` that estimates the
#' ideal `ribosome_penalty` value for the input data.
#'
#' Based on observations of true single cell data, we find that ribosomal RNA
#' loss occurs less frequently than expected based on abundance alone. To
#' adjust for this, the probability scores of ribosomal gene loss are multiplied
#' by a numerical value (`ribosome_penalty`) between 0 and 1. Lower values
#' (closer to zero) better approximate true data, with a default of 0.01,
#' though this can often be greatly refined for the input data.
#'
#' Refinement follows a similar workflow to `detect_damage`, but rather than
#' evaluating the similarity of true cells to sets of artificial cells to
#' infer their level of damage, we evaluate the similarity of artificial cells
#' to true cells to infer the effectiveness of their approximation to true
#' data. This is calculated using the distance to the nearest true cell (dTNN)
#' taken for each artificial cell found using the Euclidean distance matrix.
#' The median dTNN is computed iteratively until stabilization or a worsening
#' trend. The ideal `ribosomal_penalty` is then selected as that which
#' generated the lowest dTNN.
#'
#' @inheritParams simulate_counts
#' @param mito_quantile Numeric specifying below what proportion of
#'  mitochondrial content cells are used for sampling for simulation.
#'
#' * Default is 0.75, meaning only cells with less than 0.75 proportion of
#' mitochondrial counts are sampled for simulated.
#' @param penalty_range Numerical vector of length 2 specifying the lower
#'  and upper limit of values tested for the ribosomal penalty.
#'
#'  * Default is c(0.00001, 0.5).
#' @param penalty_step Numeric specifying the value added to each increment
#'  of penalty tested.
#'
#'  * Default is 0.005.
#' @param max_penalty_trials Numeric specifying the maximum number of
#'   iterations for the ribosomal penalty value.
#'
#'   * Default is 10.
#' @param stability_limit Numeric specifying the number of additional iterations
#'  allotted after the median minimum distance of the artificial cells to the
#'  true cells is greater than the previous minimum distance.
#'
#'  The idea here is that if a higher penalty is not causing an improvement
#'  in the output, there is little need to continue testing with larger
#'  penalties.
#'
#'  * Default is 3.
#' @param return_output String specifying what form the output of the function
#'  should take where the options are either,
#'
#'  * "penalty"
#'  * "full"
#'
#'  "Penalty" will return only the ribosomal penalty that resulted in the
#'  best performance (the smallest median distance between artificial and
#'  true cells). While "full" will return the ideal ribosomal penalty and
#'  the median distance between artificial and true cells for each penalty
#'  tested. This allows insight into how the penalty was selected.
#'
#' * Default is "penalty".
#' @param verbose Boolean specifying whether messages and function progress
#'  should be displayed in the console.
#'
#'  * Default is TRUE.
#' @return Numeric representing the ideal ribosomal penalty for an input dataset.
#' @importFrom dplyr %>% mutate if_else
#' @importFrom stats quantile
#' @importFrom scales rescale
#' @export
#' @examples
#' data("test_counts", package = "DamageDetective")
#'
#' penalty <- select_penalty(
#'  count_matrix = test_counts,
#'  stability_limit = 1,
#'  max_penalty_trials = 1,
#'  seed = 7
#' )
select_penalty <- function(
    count_matrix,
    organism = "Hsap",
    mito_quantile = 0.75,
    penalty_range = c(0.00001, 0.5),
    penalty_step = 0.005,
    max_penalty_trials = 10,
    target_damage = c(0.2, 0.99),
    damage_distribution = "right_skewed",
    distribution_steepness = "steep",
    beta_shape_parameters = NULL,
    stability_limit = 3,
    damage_proportion = 0.15,
    annotated_celltypes = FALSE,
    return_output = "penalty",
    ribosome_penalty = NULL,
    seed = NULL,
    verbose = TRUE
) {
  # Data preparation
  .check_penalty_inputs(
    mito_quantile, penalty_range, penalty_step,
    max_penalty_trials, stability_limit, return_output
  )

  gene_idx <- get_organism_indices(count_matrix, organism)
  df <- .prepare_penalty_data(count_matrix, gene_idx, mito_quantile)

  # Subset low mitochondrial proportion cells
  filtered_cells <- rownames(df)[df$mito <= quantile(df$mito, probs = mito_quantile, na.rm = TRUE)]
  filtered_matrix <- count_matrix[, filtered_cells]

  # Iterate through penalties to generate artificial cells
  penalties <- seq(penalty_range[[1]], penalty_range[[2]], by = penalty_step)
  penalty_results <- .iterate_penalties(
    penalties, filtered_matrix, df, target_damage, damage_distribution,
    distribution_steepness, beta_shape_parameters, damage_proportion,
    stability_limit, max_penalty_trials, seed, verbose
  )

  # Return the selected penalty or full results
  return(.finalise_penalty_output(penalty_results, return_output))
}

.check_penalty_inputs <- function(
    mito_quantile, penalty_range, penalty_step,
    max_penalty_trials, stability_limit, return_output
) {
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

.prepare_penalty_data <- function(count_matrix, gene_idx, mito_quantile) {
  mito <- colSums(count_matrix[gene_idx$mito_idx, , drop = FALSE]) / colSums(count_matrix)
  ribo <- colSums(count_matrix[gene_idx$ribo_idx, , drop = FALSE]) / colSums(count_matrix)
  features <- colSums(count_matrix != 0)

  df <- data.frame(
    mito = mito,
    ribo = ribo,
    features = features,
    row.names = colnames(count_matrix)
  )

  mito_top <- quantile(df$mito, probs = mito_quantile, na.rm = TRUE)
  df$mito_label <- ifelse(df$mito > mito_top, "high_mito", "low_mito")
  df$damage_label <- "cell"
  df$Damaged_Level <- 0

  return(df)
}

.iterate_penalties <- function(
    penalties, filtered_matrix, df, target_damage, damage_distribution,
    distribution_steepness, beta_shape_parameters, damage_proportion,
    stability_limit, max_penalty_trials, seed, verbose
) {
  penalty_results <- data.frame(
    Penalty = numeric(0),
    Global_mean = numeric(0)
  )

  best_dTNN_mean <- Inf
  stability_counter <- 0
  penalty_count <- 0

  for (penalty in penalties) {
    if (penalty_count >= max_penalty_trials) {
      message("Maximum penalty trials reached (", max_penalty_trials, "). Stopping.")
      break
    }

    penalty_count <- penalty_count + 1

    if (verbose) {
      message("Testing penalty of ", penalty, "...")
    }

    # Run simulation
    penalty_output <- simulate_counts(
      count_matrix = filtered_matrix,
      damage_proportion = damage_proportion,
      target_damage = target_damage,
      ribosome_penalty = penalty,
      damage_distribution = damage_distribution,
      distribution_steepness = distribution_steepness,
      seed = seed,
      generate_plot = FALSE
    )

    # Extract QC metrics and compute dTNN
    df_damaged <- .extract_damaged_cells(penalty_output$qc_summary)
    combined_df <- .combine_true_and_artificial_cells(df, df_damaged)
    dTNN_means <- .compute_dTNN(combined_df)

    # Store results
    penalty_results <- .store_penalty_results(
      penalty_results, penalty, dTNN_means, stability_counter,
      best_dTNN_mean, stability_limit
    )

    # Update best dTNN mean and stability counter
    best_dTNN_mean <- min(best_dTNN_mean, mean(dTNN_means, na.rm = TRUE))
    if (mean(dTNN_means, na.rm = TRUE) >= best_dTNN_mean) {
      stability_counter <- stability_counter + 1
      if (stability_counter > stability_limit) {
        message("Stopping early: dTNN is no longer improving.")
        break
      }
    } else {
      stability_counter <- 0
    }
  }

  return(penalty_results)
}

.extract_damaged_cells <- function(qc_summary) {
  df_damaged <- qc_summary %>%
    dplyr::filter(.data$Damaged_Level != 0)

  colnames(df_damaged)[colnames(df_damaged) == "New_RiboProp"] <- "ribo"
  colnames(df_damaged)[colnames(df_damaged) == "New_MitoProp"] <- "mito"
  colnames(df_damaged)[colnames(df_damaged) == "New_Features"] <- "features"

  df_damaged <- df_damaged[, c("ribo", "mito", "features", "Damaged_Level")]
  df_damaged$damage_label <- "artificial"
  df_damaged$mito_label <- "-"

  return(df_damaged)
}

.combine_true_and_artificial_cells <- function(df, df_damaged) {
  combined_df <- rbind(df, df_damaged)

  pca_result <- prcomp(
    combined_df[, c("mito", "ribo", "features")],
    center = TRUE, scale. = TRUE
  )
  combined_df$pca1 <- pca_result$x[, 1]
  combined_df$pca2 <- pca_result$x[, 2]

  return(combined_df)
}

.compute_dTNN <- function(combined_df) {
  true_cells <- combined_df %>%
    dplyr::filter(.data$damage_label == "cell")
  artificial_cells <- combined_df %>%
    dplyr::filter(.data$damage_label == "artificial")

  dists <- sqrt((outer(artificial_cells$pca1, true_cells$pca1, "-")^2) +
                  (outer(artificial_cells$pca2, true_cells$pca2, "-")^2))
  nearest_dists <- apply(dists, 1, min)

  artificial_cells$nearest_distance <- nearest_dists
  artificial_cells$damage_level <- cut(
    artificial_cells$Damaged_Level,
    breaks = seq(0, 1, by = 0.1),
    include.lowest = TRUE,
    labels = FALSE
  )

  dTNN_means <- tapply(
    artificial_cells$nearest_distance, artificial_cells$damage_level, mean
  )

  return(dTNN_means)
}

.store_penalty_results <- function(
    penalty_results, penalty, dTNN_means, stability_counter,
    best_dTNN_mean, stability_limit
) {
  dTNN_means_df <- as.data.frame(t(dTNN_means))
  colnames(dTNN_means_df) <- paste0("Mean_", names(dTNN_means))

  new_row <- data.frame(
    Penalty = penalty,
    Global_mean = mean(dTNN_means, na.rm = TRUE),
    dTNN_means_df
  )

  penalty_results <- dplyr::bind_rows(penalty_results, new_row)
  return(penalty_results)
}

.finalise_penalty_output <- function(penalty_results, return_output) {
  best_penalty <- penalty_results$Penalty[which.min(penalty_results$Global_mean)]

  if (return_output == "penalty") {
    return(best_penalty)
  }

  if (return_output == "full") {
    return(list(
      penalty_results = penalty_results,
      selected_penalty = best_penalty
    ))
  }
}
