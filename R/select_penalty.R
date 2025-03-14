#' select_penalty
#'
#' @param count_matrix Matrix or dgCMatrix containing the counts from
#'  single cell RNA sequencing data.
#' @param organism  String specifying the organism of origin of the input
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
#'  tested. This allows insight into how the penalty was seleected.
#'
#' * Default is "penalty".
#' @param generate_plot Boolean specifying whether the quality control plots
#'  of each penalty tested should be stored for viewing.
#'
#' * Default is FALSE.
#' @param verbose Boolean specifying whether messages and function progress
#'  should be displayed in the console.
#'
#'  * Default is TRUE.
#' @return Numeric representing the ideal ribosomal penalty for an input dataset.
#' @export
#' @examples
#' data("test_counts", package = "DamageDetective")
#'
#' penalty <- select_penalty(
#' matrix = test_counts
#' )
#' @importFrom dplyr %>% mutate if_else
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn labs theme_classic theme element_rect element_text
#' @importFrom scales rescale
select_penalty <- function(
    count_matrix,
    organism = "Hsap",
    mito_quantile = 0.75,
    penalty_range = c(0.001, 0.5),
    penalty_step = 0.005,
    stability_limit = 3,
    return_output = "penalty",
    generate_plot = FALSE,
    verbose = TRUE,
    ...
)
{
  # Phase One: Data preparations ----

  # Retrieve genes corresponding to the organism of interest

  # Human
  if (organism == "Hsap") {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^(RPS|RPL)"
  }

  # Mouse
  if (organism == "Mmus") {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^(rps|rpl)"
  }

  # Allow for user specification for non-standard organism
  if (!organism %in% c("Hsap", "Mmus")) {
    mito_pattern <- organism$mito_pattern
    ribo_pattern <- organism$ribo_pattern
  }

  # Isolate gene set indices (consistent across cells, not subsetting the matrix)
  mito_idx <- grep(mito_pattern, rownames(count_matrix), ignore.case = FALSE)
  ribo_idx <- grep(ribo_pattern, rownames(count_matrix), ignore.case = FALSE)
  mito <- colSums(count_matrix[mito_idx, , drop = FALSE]) / colSums(count_matrix)
  ribo <- colSums(count_matrix[ribo_idx, , drop = FALSE]) / colSums(count_matrix)
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
      damage_proportion = 0.15,
      target_damage = c(0.1, 0.99),
      ribosome_penalty = penalty,
      damage_distribution = "symmetric",
      damage_steepness = "steep",
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


    # Phase Three: PCA to compare penalty iterations ----
    # Best iteration must minimise the PC distance between artificial cells and true cells

    # Run PCA between all artificial and true cells
    pca_result <- prcomp(combined_df[, c("mito", "ribo", "features")], center = TRUE, scale. = TRUE)
    combined_df$pca1 <- pca_result$x[, 1]
    combined_df$pca2 <- pca_result$x[, 2]

    # Separate true and artificial cells for distance calculation
    true_cells <- subset(combined_df, damage_label == "cell")
    artificial_cells <- subset(combined_df, damage_label == "artificial")

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
      plot <- ggplot(combined_df, aes(x = ribo, y = mito)) +
        geom_point(aes(color = nearest_distance)) +
        scale_color_gradientn(
          colours = c("lightgray", "#7023FD", "#E60006"),
          values = scales::rescale(c(0, 0.1, 1)),
          limits = c(0, max(combined_df$nearest_distance, na.rm = TRUE))
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
    return(selected_penalty = best_penalty)
  }

  if (return_output == "full"){
    return(list(penalty_results = penalty_results,
                selected_penalty = best_penalty))
  }

}

