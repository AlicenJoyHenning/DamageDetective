#' detect_damage
#'
#' Quality control function to identify and filter damaged cells from an
#' input count matrix, where 'damage' is defined by the loss of cytoplasmic RNA.
#'
#' Using the simulation framework of `simulate_counts()`, `detect_damage()`
#' generates artificially damaged cell profiles by introducing defined levels
#' of RNA loss into the input data. True and artificial cells are then
#' merged and pre-processed to compute the following quality control metrics:
#'
#' * Log-normalized feature count
#' * Log-normalized total counts
#' * Mitochondrial proportion
#' * Ribosomal proportion
#' * Log-normalized MALAT1 gene expression
#'
#' Principal component analysis (PCA) is performed on these metrics,
#' and a Euclidean distance matrix is constructed from the PC embeddings.
#' For each true cell, the proportion of nearest neighbors that are
#' artificial cells (pANN) is calculated across all damage levels and the
#' damage level with the highest pANN is assigned to the true cell.
#' Finally, cells exceeding a specified damage threshold, `filter_threshold`,
#' are marked as damaged.
#'
#' This filtering method is inspired by approaches developed for DoubletFinder
#' (McGinnis et al., 2019) to detect doublets in single-cell data.
#'
#' @references
#' McGinnis, C. S., Murrow, L. M., & Gartner, Z. J. (2019). DoubletFinder:
#' Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial
#' Nearest Neighbors. *Cell Systems, 8*(4), 329-337.e4. \doi{10.1016/j.cels.2019.03.003}
#'
#' @inheritParams simulate_counts
#' @param filter_threshold Numeric specifying the proportion of RNA loss
#'  above which a cell should be considered damaged.
#'
#'  * Default is 0.75.
#' @param project_name String specifying a project identifier. Intended
#'  for generating quality control plots that are distinct across samples.
#'  However, is otherwise not relevant to the function.
#'
#'  * Default is "Project"
#' @param damage_levels Numeric specifying the number of distinct sets of
#'  artificial damaged cells simulated, each with a defined range of loss.
#'  Default ptions include,
#'
#'  * 3 : c(0.00001, 0.08), c(0.1, 0.4), c(0.5, 0.9)
#'  * 5 : c(0.00001, 0.08), c(0.1, 0.3), c(0.3, 0.5), c(0.5, 0.7), c(0.7, 0.9)
#'  * 7 : c(0.00001, 0.08), c(0.1, 0.3), c(0.3, 0.4), c(0.4, 0.5), c(0.5, 0.7),
#'  c(0.7, 0.9), c(0.9, 0.99999).
#'
#'  A user can also provide a list specifying sets with their own ranges of loss,
#'
#'  * damage_levels = list(
#'    pANN_50 = c(0.1, 0.5),
#'    pANN_100 = c(0.5, 1)
#'  )
#'
#'  By introducing more sets of damage a user can improve the accuracy of
#'  loss estimations (scaled_pANN) as they are found through scaling the pANN
#'  within each set according to the lower and upper boundary of the set's
#'  damage level. However, introducing more sets increases the computational
#'  time for the function.
#'
#' * Default is 5.
#' @param mito_quantile Numeric between 0 and 1 specifying below what
#'  level of mitochondrial proportion cells are sampled for simulations.
#'  This step is done to protect against simulating damaged cell profiles
#'  from cells that are likely damaged.
#'
#'  * Default is 0.75.
#' @param kN Numeric describing how many nearest neighbours are considered
#'  for pANN calculations. kN cannot exceed the total cell number.
#'
#'  * Default is one third of the total cell number.
#' @param filter_counts Boolean specifying whether the output matrix
#'  should be filtered, returned containing only cells that fall below
#'  the filter threshold. Alternatively, a data frame containing cell
#'  barcodes and their associated label as either 'damaged' or 'cell'
#'  is returned.
#'
#'  * Default is FALSE.
#' @param verbose Boolean specifying whether messages and function progress
#'  should be displayed in the console.
#'
#'  * Default is TRUE.
#' @return Filtered matrix or data frame containing damage labels.
#' @importFrom Matrix colSums Matrix
#' @importFrom RANN nn2
#' @importFrom dplyr %>% group_by summarise mutate case_when first pull
#' @importFrom ggplot2 ggsave theme element_rect margin
#' @importFrom cowplot ggdraw draw_label plot_grid
#' @importFrom ggpubr get_legend
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom Seurat ScaleData PercentageFeatureSet FetchData
#' @importFrom rlang expr !!! .data
#' @importFrom stats prcomp
#' @export
#' @examples
#' data("test_counts", package = "DamageDetective")
#'
#' test <- detect_damage(
#'   count_matrix = test_counts,
#'   ribosome_penalty = 0.001,
#'   project_name = "Test",
#'   generate_plot = FALSE
#' )
detect_damage <- function(
    count_matrix,
    ribosome_penalty,
    organism = "Hsap",
    annotated_celltypes = FALSE,
    target_damage = c(0.1, 0.8),
    damage_distribution = "right_skewed",
    distribution_steepness = "moderate",
    beta_shape_parameters = NULL,
    project_name = "Project",
    filter_threshold = 0.75,
    damage_levels = 5,
    damage_proportion = 0.15,
    mito_quantile = 0.75,
    kN = NULL,
    generate_plot = TRUE,
    filter_counts = FALSE,
    verbose = TRUE
){
  # Output warnings if user inputs are invalid ----
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

  # Data preparations for simulation ----

  # Retrieve genes corresponding to the organism of interest
  gene_idx <- get_organism_indices(count_matrix, organism)

  # Collect mito & ribo information of all true cells
  mito <- colSums(count_matrix[gene_idx$mito_idx, , drop = FALSE]) / colSums(count_matrix)
  ribo <- colSums(count_matrix[gene_idx$ribo_idx, , drop = FALSE]) / colSums(count_matrix)

  # Combine mito and ribo into a single matrix and preserve column names
  mito_ribo_data <- cbind(mito, ribo)
  colnames(mito_ribo_data) <- c("mito", "ribo")
  rownames(mito_ribo_data) <- colnames(count_matrix)

  # Identify distribution of mitochondrial proportion for simulation sampling
  mito_top <- quantile(mito_ribo_data[, "mito"], probs = mito_quantile, na.rm = TRUE)

  # Subset true cells below the cells in the top distribution using Matrix operations
  filtered_cells <- rownames(mito_ribo_data)[mito_ribo_data[, "mito"] <= mito_top]

  # Efficiently subset count_matrix while preserving sparsity
  matrix_filtered <- count_matrix[, filtered_cells, drop = FALSE]

  # Generate simulations for defined damage levels ----

  # Retrieve damage levels for simulation
  if (is.list(damage_levels)) {
    ranges <- damage_levels
  }

  if (damage_levels == 3) {
    ranges <- list(
      pANN_10 = c(0.00001, 0.08),
      pANN_40 = c(0.1, 0.4),
      pANN_90 = c(0.5, 0.9)
    )
  }

  if (damage_levels == 5) {
    ranges <- list(
      pANN_10 = c(0.00001, 0.08),
      pANN_30 = c(0.1, 0.3),
      pANN_50 = c(0.3, 0.5),
      pANN_70 = c(0.5, 0.7),
      pANN_90 = c(0.7, 0.9)
    )
  }

  if (damage_levels == 7) {
    ranges <- list(
      pANN_10 = c(0.00001, 0.08),
      pANN_30 = c(0.1, 0.3),
      pANN_40 = c(0.3, 0.4),
      pANN_50 = c(0.4, 0.5),
      pANN_70 = c(0.5, 0.7),
      pANN_90 = c(0.7, 0.9),
      pANN_100 = c(0.9, 0.99999)
    )
  }

  # Initialize list to store damaged matrices
  damaged_matrices <- list()
  damage_ranges <- c()

  # Iterate through ranges and simulate damage
  for (i in seq_along(ranges)) {
    level <- names(ranges)[i]
    damage_range <- ranges[[i]]
    damage_ranges <- append(damage_range, damage_ranges)

    if (verbose) {
      message("Simulating cells between ", damage_range[[1]], " and ", damage_range[[2]], " RNA loss...")
    }

    # Simulate counts
    damaged_cells <- simulate_counts(
      count_matrix = matrix_filtered,
      damage_proportion = damage_proportion,
      target_damage = damage_range,
      ribosome_penalty = ribosome_penalty,
      damage_distribution = damage_distribution,
      distribution_steepness = distribution_steepness,
      generate_plot = FALSE
    )

    # Extract barcodes of damaged cells
    barcodes <- damaged_cells$qc_summary %>%
      dplyr::filter(.data$Damaged_Level != 0) %>%
      dplyr::pull(.data$Cell)

    damaged_matrix <- damaged_cells$matrix[, barcodes]

    # Append the level of damage to the barcodes
    colnames(damaged_matrix) <- paste0(colnames(damaged_matrix), "_", gsub("pANN_", "", level))

    # Store the damaged matrix in the list
    damaged_matrices[[level]] <- damaged_matrix
  }

  # Combine all damaged matrices at once
  matrix_updated <- do.call(cbind, damaged_matrices)


  # Combine and remove duplicated, unaltered true cells
  matrix_combined <- cbind(matrix_updated, count_matrix)
  matrix_combined <- Matrix::Matrix(matrix_combined, sparse = TRUE)
  matrix_combined <- matrix_combined[, unique(colnames(matrix_combined))]


  # Compute quality control measures for the matrix
  features <- Matrix::colSums(matrix_combined > 0)
  log.features <- log10(features)
  total_counts <- Matrix::colSums(matrix_combined)
  log.counts <- log10(total_counts)

  # Compute mitochondrial and ribosomal proportions
  mt.prop <- Matrix::colSums(matrix_combined[grep(gene_idx$mito_pattern, rownames(matrix_combined)), , drop = FALSE]) / total_counts
  rb.prop <- Matrix::colSums(matrix_combined[grep(gene_idx$ribo_pattern, rownames(matrix_combined)), , drop = FALSE]) / total_counts

  # Avoid division by zero
  mt.prop[is.na(mt.prop)] <- 0
  rb.prop[is.na(rb.prop)] <- 0

  # Extract MALAT1 expression and log-transform
  malat1_expr <- matrix_combined[gene_idx$MALAT1, , drop = FALSE]
  malat1 <- as.vector(malat1_expr)

  # Extract Damage_level from cell names
  damage_numbers <- sub(".*_", "", names(ranges))
  cell_names <- colnames(matrix_combined)
  Damage_level <- sub(".*_", "", cell_names)
  Damage_level <- ifelse(!Damage_level %in% damage_numbers, 0, Damage_level)

  # Combine into a dataframe
  metadata_stored <- data.frame(
    features = features,
    log.features = log.features,
    counts = total_counts,
    log.counts = log.counts,
    mt.prop = mt.prop,
    rb.prop = rb.prop,
    malat1 = malat1,
    Damage_level = Damage_level,
    row.names = cell_names
  )


  # Perform PCA on the quality control statistics of all cells ----

  # Isolate variables of interest
  metadata_pca <- metadata_stored[, c("log.features", "log.counts", "mt.prop", "rb.prop", "malat1")]

  # Perform PCA
  pca_result <- prcomp(metadata_pca, center = TRUE, scale. = TRUE)

  # Extract PCA embeddings of top principal components
  pca_coord <- pca_result$x[, 1:5]

  # Perform nearest neighbor search to get kN nearest neighbors

  # If number of neighbours is not specified, default to third of total
  if (is.null(kN)){
    kN <- round(dim(count_matrix)[2] / 5, 0)
  }

  # nn2 function: k = number of neighbors, search on pca_coord
  knn_result <- RANN::nn2(pca_coord, pca_coord, k = kN)

  # Get the neighbor indices (without the cell itself)
  neighbor_indices <- knn_result$nn.idx

  # Isolate columns & ensure cell names present
  metadata_plot <- metadata_stored[, c("features", "counts", "mt.prop", "rb.prop", "malat1", "Damage_level")]
  metadata_plot$Cells <- rownames(metadata_stored)

  # Define sets of cells based on damage level
  barcodes <- list()
  for (number in damage_numbers){
    barcode_name <- paste0("barcode_", number)
    barcodes[[barcode_name]] <- metadata_plot$Cells[metadata_plot$Damage_level == number]
  }

  # Proportion of damage level-specific artificial nearest neighbours (pANNs)
  compute_pANN <- function(barcode_set) {
    sapply(rownames(metadata_stored), function(cell) {

      # Get the neighbors for the cell (excluding itself)
      neighbors <- neighbor_indices[which(rownames(metadata_stored) == cell), -1]
      neighbor_barcodes <- rownames(metadata_stored)[neighbors]

      # Compute proportion of neighbors in the barcode set
      sum(neighbor_barcodes %in% barcode_set) / kN
    })
  }

  # Compute pANN for different damage levels
  if (verbose){
    message("Computing pANN...")
  }

  # Compute proportion of artificial neighbours
  pANN_names <- c()
  for (number in damage_numbers){
    barcode_name <- paste0("barcode_", number)
    pANN_name <- paste0("pANN_", number)
    pANN_names <- append(pANN_name, pANN_names)
    metadata_plot[[pANN_name]] <- compute_pANN(barcodes[[barcode_name]])
  }

  # Find which damage population a cell is most frequently neighboring

  # Isolate true cells & define relevant columns
  metadata_plot <- metadata_plot %>%
    dplyr::filter(.data$Damage_level == 0)

  # Find the column name with the maximum value for each true cell
  metadata_plot$max_pANN_col <- apply(metadata_plot[pANN_names], 1, function(row) {
    pANN_names[which.max(row)]
  })

  # Calculate min/max pANN value for each damage level
  pANN_summary_stats <- metadata_plot %>%
    dplyr::group_by(.data$max_pANN_col) %>%
    dplyr::summarise(
      min_value = min(get(first(.data$max_pANN_col)), na.rm = TRUE),
      max_value = max(get(first(.data$max_pANN_col)), na.rm = TRUE)
    )

  # Use the min/max values to scale pANN of each cell to lie between range
  metadata_plot <- metadata_plot %>%
    mutate(
      lower_scale = case_when(
        !!!lapply(names(ranges), function(name) {
          expr(.data$max_pANN_col == !!name ~ ranges[[!!name]][[1]])
        })
      ),
      upper_scale = case_when(
        !!!lapply(names(ranges), function(name) {
          expr(.data$max_pANN_col == !!name ~ ranges[[!!name]][[2]])
        })
      ),
      min = case_when(
        !!!lapply(names(ranges), function(name) {
          expr(.data$max_pANN_col == !!name ~ pANN_summary_stats$min_value[pANN_summary_stats$max_pANN_col == !!name])
        })
      ),
      max = case_when(
        !!!lapply(names(ranges), function(name) {
          expr(.data$max_pANN_col == !!name ~ pANN_summary_stats$max_value[pANN_summary_stats$max_pANN_col == !!name])
        })
      ),
      value = case_when(
        !!!lapply(names(ranges), function(name) {
          expr(.data$max_pANN_col == !!name ~ get(!!name))
        })
      )
    )

  # Apply scaling
  metadata_plot$DamageDetective <- metadata_plot$lower_scale + ((metadata_plot$value - metadata_plot$min) / (metadata_plot$max - metadata_plot$min)) * (metadata_plot$upper_scale - metadata_plot$lower_scale)

  # Identify and visualize damaged cells ----

  metadata_plot$DamageDetective_filter <- ifelse(metadata_plot$DamageDetective >= filter_threshold, "damaged", "cell")
  metadata_output <- metadata_plot[, c("Cells", "DamageDetective")]

  # If specified, filter and return the count matrix only
  if (filter_counts){
    metadata_filtered <- metadata_plot %>%
      dplyr::filter(.data$DamageDetective_filter == "cell")
    final_filtered_cells <- metadata_filtered$Cells
    final_filtered_matrix <- count_matrix[, final_filtered_cells]
    output <- final_filtered_matrix
  } else {
    output <- metadata_output
  }

  # Visualise cells according to damage level
  if (generate_plot){

    # Generate QC plots with cells coloured by scaled pANN
    plot_mito_ribo <- plot_outcome(
      data = metadata_plot,
      x = "rb.prop",
      y = "mt.prop",
      damage_column = "DamageDetective",
      altered = TRUE,
      mito_ribo = TRUE,
      target_damage = c(0.5, 0.9)
    )

    plot_mito_features <- plot_outcome(
      data = metadata_plot,
      x = "features",
      y = "mt.prop",
      damage_column = "DamageDetective",
      altered = TRUE,
      target_damage = c(0.5, 0.9)
    )

    # Extract the legend from mito_ribo_new
    legend <- ggpubr::get_legend(plot_mito_features)

    # Create title for the plot
    title_text <- paste0(project_name, " damage level predictions")
    title <- cowplot::ggdraw() +
      cowplot::draw_label(title_text, fontface = 'bold', hjust = 0.5)

    # Arrange altered plots in a single row
    mito_ribo_no_legend <- plot_mito_ribo + ggplot2::theme(legend.position = "none")
    mito_features_no_legend <- plot_mito_features + ggplot2::theme(legend.position = "none")
    combined_plots <- cowplot::plot_grid(mito_features_no_legend, mito_ribo_no_legend, ncol = 2)
    final_plot <- cowplot::plot_grid(title, combined_plots, legend, ncol = 1, rel_heights = c(0.2, 1, 0.25))

    # Increase margins around the total plot area
    final_plot <- final_plot +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5, 20, 20, 20),
        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
        plot.background = ggplot2::element_rect(fill = "white", color = "white")
      )

    # Display the final plot
    print(final_plot)

    output <- list(
      output = output,
      plot = final_plot
    )
  }

  return(output)
}
