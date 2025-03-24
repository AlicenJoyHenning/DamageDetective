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
#' @importFrom fields rdist
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
  # Ensure user inputs are valid ----

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

  # Ensure count matrix is of 'matrix' form
  count_matrix <- as.matrix(count_matrix)

  # Data preparations for simulation ----

  # Retrieve genes corresponding to the organism of interest
  gene_idx <- get_organism_indices(count_matrix, organism)

  # Collect mito & ribo information of all true cells
  mito <- colSums(count_matrix[gene_idx$mito_idx, , drop = FALSE]) / colSums(count_matrix)
  ribo <- colSums(count_matrix[gene_idx$ribo_idx, , drop = FALSE]) / colSums(count_matrix)

  df <- data.frame(mito = mito,
                   ribo = ribo,
                   row.names = colnames(count_matrix))

  # Identify high mito proportion cells (likely damaged)
  mito_top <- quantile(df$mito, probs = mito_quantile, na.rm = TRUE)

  # Subset true cells below the 75th percentile
  filtered_cells <- rownames(df)[df$mito <= mito_top]
  matrix_filtered <- count_matrix[, filtered_cells]

  # Generate simulations for defined damage levels ----

  # Retrieve damage levels for simulation
  if (is.list(damage_levels)){
    ranges <- damage_levels
  }

  if (damage_levels == 3){
    ranges = list(
      pANN_10 = c(0.00001, 0.08),
      pANN_40 = c(0.1, 0.4),
      pANN_90 = c(0.5, 0.9)
    )
  }

  if (damage_levels == 5){
    ranges = list(
      pANN_10 = c(0.00001, 0.08),
      pANN_30 = c(0.1, 0.3),
      pANN_50 = c(0.3, 0.5),
      pANN_70 = c(0.5, 0.7),
      pANN_90 = c(0.7, 0.9)
    )
  }

  if (damage_levels == 7){
    ranges = list(
      pANN_10 = c(0.00001, 0.08),
      pANN_30 = c(0.1, 0.3),
      pANN_40 = c(0.3, 0.4),
      pANN_50 = c(0.4, 0.5),
      pANN_70 = c(0.5, 0.7),
      pANN_90 = c(0.7, 0.9),
      pANN_100 = c(0.9, 0.99999)
    )
  }

  # Initialize matrix with original counts
  matrix_updated <- matrix_filtered

  # Iterate through ranges and simulate damage
  damage_ranges <- c()
  for (i in seq_along(ranges)) {
    level <- names(ranges)[i]
    damage_range <- ranges[[i]]
    damage_ranges <- append(ranges[[i]], damage_ranges)

    if (verbose){
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
      generate_plot = generate_plot
    )

    # Extract barcodes of damaged cells
    barcodes <- damaged_cells$qc_summary %>%
      dplyr::filter(.data$Damaged_Level != 0) %>%
      dplyr::pull(.data$Cell)

    damaged_matrix <- damaged_cells$matrix[, barcodes]

    # Append the level of damage to the barcodes
    colnames(damaged_matrix) <- paste0(colnames(damaged_matrix), "_", gsub("pANN_", "", level))

    # Merge with updated matrix
    matrix_updated <- cbind(matrix_updated, damaged_matrix)

  }

  # Combine and remove duplicated, unaltered true cells
  matrix_combined <- cbind(matrix_updated, count_matrix)
  matrix_combined <- matrix_combined[, unique(colnames(matrix_combined))]

  # Create Seurat object with all cells
  combined <- suppressWarnings(CreateSeuratObject(counts = matrix_combined))

  # Pre-process the Seurat object
  combined <- NormalizeData(combined, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE)

  # Add meta_data label with level of damage of the cell
  damage_numbers <- sub(".*_", "", names(ranges))
  combined$Damage_level <- sub(".*_", "", rownames(combined@meta.data))
  combined$Damage_level <- ifelse(!combined$Damage_level %in% damage_numbers, 0, combined$Damage_level)

  # Find mitochondrial and ribosomal percentage for each cell
  combined$log.features <- log10(combined$nFeature_RNA)
  combined$log.counts <- log10(combined$nCount_RNA)
  combined$mt.prop <- PercentageFeatureSet(combined, pattern = gene_idx$mito_pattern) / 100
  combined$rb.prop <- PercentageFeatureSet(combined, pattern = gene_idx$ribo_pattern) / 100

  # Avoid infinite values (log10(0))
  combined$malat1 <- FetchData(combined, vars = gene_idx$MALAT1)
  epsilon <- 1e-6
  combined$malat1 <- log10(combined$malat1 + epsilon)

  # Perform PCA on the quality control statistics of all cells ----

  # Isolate variables of interest
  metadata <- combined@meta.data
  metadata <- metadata[, c("log.features", "log.counts", "mt.prop", "rb.prop", "malat1")]

  # Perform PCA
  pca_result <- prcomp(metadata, center = TRUE, scale. = TRUE)

  # Extract PCA embeddings of top principal components
  pca_coord <- pca_result$x[, 1:5]

  # Calculate the euclidean distance between PC embeddings of cells
  dist_mat <- rdist(pca_coord)  # Euclidean distances
  dimnames(dist_mat) <- list(rownames(metadata), rownames(metadata))

  # Isolate columns of for PCA & ensure cell names present
  metadata <- combined@meta.data[, c("nFeature_RNA", "nCount_RNA", "mt.prop", "rb.prop", "malat1", "Damage_level")]
  metadata$Cells <- rownames(metadata)

  # Define sets of cells based on damage level
  barcodes <- list()
  for (number in damage_numbers){
    barcode_name <- paste0("barcode_", number)
    barcodes[[barcode_name]] <- metadata$Cells[metadata$Damage_level == number]
  }

  # If number of neighbours is not specified, default to third of total
  if (is.null(kN)){
    kN <- round(dim(count_matrix)[2] / 3, 0)
  }

  # Proportion of damage level-specific artificial nearest neighbours (pANNs)
   compute_pANN <- function(barcode_set) {
    sapply(rownames(dist_mat), function(cell) {
      # Find indices of the nearest neighbors (excluding itself)
      kNN <- kN + 1
      neighbors <- order(dist_mat[cell, ])[2:kN]
      neighbor_barcodes <- rownames(dist_mat)[neighbors]

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
    metadata[[pANN_name]] <- compute_pANN(barcodes[[barcode_name]])
  }

  # Find which damage population a cell is most frequently neighboring

  # Isolate  true cells & define relevant columns
  metadata <- metadata %>%
    dplyr::filter(.data$Damage_level == 0)

  # Find the column name with the maximum value for each true cell
  metadata$max_pANN_col <- apply(metadata[pANN_names], 1, function(row) {
    pANN_names[which.max(row)]
  })

  # Calculate min/max pANN value for each damage level
  pANN_summary_stats <- metadata %>%
    dplyr::group_by(.data$max_pANN_col) %>%
    dplyr::summarise(
      min_value = min(get(first(.data$max_pANN_col)), na.rm = TRUE),
      max_value = max(get(first(.data$max_pANN_col)), na.rm = TRUE)
    )

  # Use the min/max values to scale pANN of each cell to lie between range
  metadata <- metadata %>%
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
  metadata$DamageDetective <- metadata$lower_scale + ((metadata$value  - metadata$min) / (metadata$max - metadata$min)) * (metadata$upper_scale - metadata$lower_scale)

  # Identify and visualise damaged cells  ----
  metadata$DamageDetective_filter <- ifelse(metadata$scaled_pANN >= filter_threshold, "damaged", "cell")
  metadata_output <- metadata[, c("Cells", "DamageDetective")]

  # If specified, filter and return the count matrix only
  if (filter_counts){
    metadata_filtered <- metadata %>%
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
      data = metadata,
      x = "rb.prop",
      y = "mt.prop",
      damage_column = "scaled_pANN",
      altered = TRUE,
      mito_ribo = TRUE,
      target_damage = c(0.5, 0.9)
    )

    plot_mito_features <- plot_outcome(
      data = metadata,
      x = "nFeature_RNA",
      y = "mt.prop",
      damage_column = "scaled_pANN",
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
