#' select_penalty
#'
#' @param matrix
#' @param mito_quantile
#' @param penalty_range
#' @param penalty_step
#' @param damage_proportion
#' @param target_damage
#' @param stability_limit
#' @param damage_distribution
#' @param return_output
#' @param damage_steepness
#' @param generate_plot
#' @param verbose
#'
#' @return numeric the value of the ideal ribosomal penalty
#' @export
#'
#' @examples
#' data("test_counts", package = "DamageDetective")
#'
#' penalty <- select_penalty(
#' matrix = test_counts
#' )
select_penalty <- function(
    matrix,
    mito_quantile = 0.75,
    penalty_range = c(0.001, 0.5),
    penalty_step = 0.005,
    damage_proportion = 0.15,
    target_damage = c(0.1, 0.99),
    stability_limit = 3,
    damage_distribution = "symmetric",
    return_output = "penalty",
    damage_steepness = "steep",
    generate_plot = FALSE,
    verbose = TRUE)
{

  # 1. Create df with mito, ribo information of all true cells
  mito_idx <- grep("^MT-", rownames(matrix), ignore.case = FALSE)
  ribo_idx <- grep("^(RPS|RPL)", rownames(matrix), ignore.case = FALSE)
  mito <- colSums(matrix[mito_idx, , drop = FALSE]) / colSums(matrix)
  ribo <- colSums(matrix[ribo_idx, , drop = FALSE]) / colSums(matrix)
  features <- colSums(matrix != 0)

  df <- data.frame(mito = mito,
                   ribo = ribo,
                   features = features,
                   row.names = colnames(matrix))

  # Identify high mito proportion cells (likely damaged)
  mito_top <- quantile(df$mito, probs = mito_quantile, na.rm = TRUE)
  df$mito_label <- ifelse(df$mito > mito_top, "high_mito", "low_mito")
  df$damage_label <- "cell"
  df$Damaged_Level <- 0

  # 2. Subset to true cells below the 75th percentile
  filtered_cells <- rownames(df)[df$mito <= mito_top]
  filtered_matrix <- matrix[, filtered_cells]

  # Define penalty values to test (starting from 0.001, increasing by 0.05, max 0.5)
  penalties <- seq(penalty_range[[1]], penalty_range[[2]], by = penalty_step)
  penalty_results <- list()
  best_dTNN_median <- Inf  # Initialize best dTNN median to a very high value

  # Define stopping criteria
  # improvement_threshold <- 0.01  # Minimum improvement needed
  stability_counter <- 0  # Counter for consecutive non-improving iterations
  max_penalty_trials <- 10  # Reduce maximum penalty values to test
  penalty_count <- 0  # Track number of trials

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
      save_plot = NULL
    )

    # Extract QC metrics from output
    df_damaged <- penalty_output$qc_summary
    df_damaged <- subset(df_damaged, Damaged_Level != 0)

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

    # 3. Perform PCA
    pca_result <- prcomp(combined_df[, c("mito", "ribo", "features")], center = TRUE, scale. = TRUE)
    combined_df$pca1 <- pca_result$x[, 1]
    combined_df$pca2 <- pca_result$x[, 2]

    # Separate true and artificial cells for distance calculation
    true_cells <- subset(combined_df, damage_label == "cell")
    artificial_cells <- subset(combined_df, damage_label == "artificial")

    # 4. Compute the nearest true cell for each artificial cell
    dists <- sqrt((outer(artificial_cells$pca1, true_cells$pca1, "-")^2) +
                    (outer(artificial_cells$pca2, true_cells$pca2, "-")^2))
    nearest_dists <- apply(dists, 1, min)  # Distance to nearest true cell

    # 5. Assign damage levels for categorization
    artificial_cells$damage_level <- cut(artificial_cells$Damaged_Level,
                                         breaks = seq(0, 1, by = 0.1),
                                         include.lowest = TRUE,
                                         labels = FALSE)
    artificial_cells$nearest_distance <- nearest_dists

    # 6. Calculate the median distance from true nearest neighbours (dTNN) by damage level
    dTNN_medians <- tapply(artificial_cells$nearest_distance, artificial_cells$damage_level, median)
    current_median_dTNN <- median(dTNN_medians)

    # Match output to cell in combined df
    combined_df$nearest_distance <- 0
    combined_df$nearest_distance[rownames(combined_df) %in% rownames(artificial_cells)] <-
      artificial_cells$nearest_distance[match(rownames(combined_df)[rownames(combined_df) %in% rownames(artificial_cells)], rownames(artificial_cells))]

    # 7. Plot the mito vs ribo for true and artificial cells, colored by dTNN
    if (generate_plot){
      plot <- ggplot(combined_df, aes(x = ribo, y = mito)) +
        geom_point(aes(color = nearest_distance)) +
        scale_color_gradientn(
          colours = c("lightgray", "#7023FD", "#E60006"),  # Light gray to purple to red
          values = scales::rescale(c(0, 0.1, 1)),  # Adjust the scale from 0 (light gray) to 1 (red)
          limits = c(0, max(combined_df$nearest_distance, na.rm = TRUE))  # Set the max distance as the upper limit
        ) +
        labs(title = paste("Penalty:", penalty, ", Median dTNN:", round(current_median_dTNN, 4)),
             y = "Mito. proportion", x = "Ribo. Prop",
             color = "dTNN") +
        theme_classic() +
        theme(panel.border = element_rect(fill = NA, colour = "black"),
              plot.title = element_text(face = "bold"))

      # Display plot
      print(plot)

      # Store results for this penalty
      penalty_results[[paste0("penalty_", penalty)]] <- list(
        median = current_median_dTNN,
        dTNN_medians = dTNN_medians,
        plot = plot
      )
    }

    # Else skip plot
    penalty_results[[paste0("penalty_", penalty)]] <- list(
      median = current_median_dTNN,
      dTNN_medians = dTNN_medians
    )

    ## 8. Stop if dTNN has stabilized
    if (length(penalty_results) == 1) {
      best_dTNN <-  current_median_dTNN
      best_penalty <- penalty
    }

    if (length(penalty_results) > 1) {
      if (current_median_dTNN < best_dTNN) {
        best_dTNN <- current_median_dTNN
        best_penalty <- penalty
        stability_counter <- 0  # Reset stability counter on improvement
      } else {
        stability_counter <- stability_counter + 1  # Increase counter when worsening

        # Stop if 2 additional trials after worsening do not improve
        if (stability_counter > stability_limit) {
          message("Stopping early: dTNN is no longer improving.")
          break
        }
      }
    }
  }

  # 9. Return the selected penalty along with penalty results
  if (return_output == "penalty"){
    return(selected_penalty = best_penalty)
  }

  if (return_output == "full"){
    return(list(penalty_results = penalty_results,
                selected_penalty = best_penalty))
  }

}
