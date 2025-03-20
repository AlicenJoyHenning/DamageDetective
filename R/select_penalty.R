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
#'  * Default is c(0.001, 0.5).
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
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn labs theme_classic theme element_rect element_text
#' @importFrom scales rescale
#' @export
#' @examples
#' data("test_counts", package = "DamageDetective")
#'
#' penalty <- select_penalty(
#'  count_matrix = test_counts,
#'  generate_plot = FALSE
#' )
select_penalty <- function(
    count_matrix,
    organism = "Hsap",
    mito_quantile = 0.75,
    penalty_range = c(0.001, 0.5),
    penalty_step = 0.005,
    max_penalty_trials = 10,
    target_damage = c(0.1, 0.99),
    damage_distribution = "right_skewed",
    damage_steepness = "steep",
    beta_shape_parameters = NULL,
    stability_limit = 3,
    damage_proportion = 0.15,
    annotated_celltypes = FALSE,
    return_output = "penalty",
    ribosome_penalty = NULL,
    generate_plot = FALSE,
    verbose = TRUE
){

  # Ensure user inputs are valid ----

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

  # Ensure count matrix is of 'matrix' form
  count_matrix <- as.matrix(count_matrix)

  # Data preparations ----

  # Retrieve genes corresponding to the organism of interest
  gene_idx <- get_organism_indices(count_matrix, organism)

  # Collect mito & ribo information
  mito <- colSums(count_matrix[gene_idx$mito_idx, , drop = FALSE]) / colSums(count_matrix)
  ribo <- colSums(count_matrix[gene_idx$ribo_idx, , drop = FALSE]) / colSums(count_matrix)
  features <- colSums(count_matrix != 0)

  df <- data.frame(mito = mito,
                   ribo = ribo,
                   features = features,
                   row.names = colnames(count_matrix))

  # Identify high mito. proportion cells (likely damaged) for simulation
  mito_top <- quantile(df$mito, probs = mito_quantile, na.rm = TRUE)
  df$mito_label <- ifelse(df$mito > mito_top, "high_mito", "low_mito")
  df$damage_label <- "cell"
  df$Damaged_Level <- 0

  # Subset low mito. proportion cells
  filtered_cells <- rownames(df)[df$mito <= mito_top]
  filtered_matrix <- count_matrix[, filtered_cells]

  # Phase Two: Iterate through penalties to generate artificial cells ----

  # Define penalty values to test (starting from 0.001, increasing by 0.05)
  penalties <- seq(penalty_range[[1]], penalty_range[[2]], by = penalty_step)
  penalty_results <- list()
  best_dTNN_median <- Inf # Initialize best dTNN median to a very high value

  # Initialize counter for consecutive non-improving iterations
  stability_counter <- 0
  max_penalty_trials <- max_penalty_trials
  penalty_count <- 0

  for (penalty in penalties) {

    if (penalty_count >= max_penalty_trials) {
      message("Maximum penalty trials reached (", max_penalty_trials, "). Stopping.")
      break
    }

    penalty_count <- penalty_count + 1  # Increase count of penalties tested

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
      damage_steepness = damage_steepness,
      generate_plot = FALSE
    )

    # Extract QC metrics from output
    df_damaged <- penalty_output$qc_summary
    df_damaged <- df_damaged %>%
      dplyr::filter(.data$Damaged_Level != 0)

    # Rename columns
    colnames(df_damaged)[colnames(df_damaged) == "New_RiboProp"] <- "ribo"
    colnames(df_damaged)[colnames(df_damaged) == "New_MitoProp"] <- "mito"
    colnames(df_damaged)[colnames(df_damaged) == "New_Features"] <- "features"

    # Rename columns and isolate the relevant columns  (mito, ribo, and damage_label)
    df_damaged <- df_damaged[, c("ribo", "mito", "features", "Damaged_Level")]

    # Assign artificial label
    df_damaged$damage_label <- "artificial"
    df_damaged$mito_label <- "-"

    # Combine true and artificial cells
    combined_df <- rbind(df, df_damaged)

    # Phase Three: PCA to compare penalty iterations ----
    # Best iteration must minimise the PC distance between artificial cells and true cells

    # Run PCA between all artificial and true cells
    pca_result <- prcomp(combined_df[, c("mito", "ribo", "features")], center = TRUE, scale. = TRUE)
    combined_df$pca1 <- pca_result$x[, 1]
    combined_df$pca2 <- pca_result$x[, 2]

    # Separate true and artificial cells for distance calculation
    true_cells <- combined_df %>%
      dplyr::filter(.data$damage_label == "cell")
    artificial_cells <- combined_df %>%
      dplyr::filter(.data$damage_label == "artificial")

    # Compute the nearest true cell for each artificial cell
    dists <- sqrt((outer(artificial_cells$pca1, true_cells$pca1, "-")^2) +
                    (outer(artificial_cells$pca2, true_cells$pca2, "-")^2))
    nearest_dists <- apply(dists, 1, min)  # Distance to nearest true cell

    # Group artificial cells into their corresponding damage level
    artificial_cells$damage_level <- cut(artificial_cells$Damaged_Level,
                                         breaks = seq(0, 1, by = 0.1),
                                         include.lowest = TRUE,
                                         labels = FALSE)
    artificial_cells$nearest_distance <- nearest_dists

    # Median PC distance from true nearest neighbours (dTNN) for each damage level
    dTNN_medians <- tapply(artificial_cells$nearest_distance, artificial_cells$damage_level, median)
    current_median_dTNN <- median(dTNN_medians)

    # Match output to cell in combined df
    combined_df$nearest_distance <- 0
    combined_df$nearest_distance[rownames(combined_df) %in% rownames(artificial_cells)] <-
      artificial_cells$nearest_distance[match(rownames(combined_df)[rownames(combined_df) %in% rownames(artificial_cells)], rownames(artificial_cells))]

    # Plot the mito vs ribo for true and artificial cells, colored by dTNN
    if (generate_plot){

      QC_plot <- plot_outcome(
        data = combined_df,
        x = "ribo",
        y = "mito",
        damage_column = "nearest_distance",
        altered = TRUE,
        mito_ribo = TRUE)

      # Scale colour by nearest distances
      plot <- QC_plot +
        scale_color_gradientn(
          colours = c("lightgray", "#7023FD", "#E60006"),
          limits = c(0, max(combined_df$nearest_distance, na.rm = TRUE))
        ) +
        labs(title = paste("Penalty:", penalty, ", Median dTNN:", round(current_median_dTNN, 4)),
             y = "Mito. proportion", x = "Ribo. Prop",
             color = "dTNN")

      # Display plot
      print(plot)

      # Store results for this penalty
      penalty_results[[paste0("penalty_", penalty)]] <- list(
        median = current_median_dTNN,
        dTNN_medians = dTNN_medians,
        plot = plot
      )
    } else {
      # Else skip plot
      penalty_results[[paste0("penalty_", penalty)]] <- list(
        median = current_median_dTNN,
        dTNN_medians = dTNN_medians
      )
    }

    # Phase Four: Stop if dTNN has stabilized ----

    if (length(penalty_results) == 1) {
      best_dTNN <-  current_median_dTNN
      best_penalty <- penalty
    }

    if (length(penalty_results) > 1) {
      if (current_median_dTNN < best_dTNN) {
        best_dTNN <- current_median_dTNN
        best_penalty <- penalty

        # Reset stability counter on improvement
        stability_counter <- 0

      } else {
        # Increase counter when worsening
        stability_counter <- stability_counter + 1

        # Stop if additional trials after worsening are not improving
        if (stability_counter > stability_limit) {
          message("Stopping early: dTNN is no longer improving.")
          break
        }
      }
    }
  }

  # Return the selected penalty along with penalty results
  if (return_output == "penalty"){
    return(best_penalty)
  }

  if (return_output == "full"){
    return(list(penalty_results = penalty_results,
                selected_penalty = best_penalty))
  }
}
