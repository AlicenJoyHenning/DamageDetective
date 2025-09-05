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
#' For each true cell, the proportion of nearest neighbours that are
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
#' Nearest neighbours. *Cell Systems, 8*(4), 329-337.e4.
#' \doi{10.1016/j.cels.2019.03.003}
#'
#' @inheritParams simulate_counts
#' @param filter_threshold Numeric specifying the proportion of RNA loss
#'  above which a cell should be considered damaged.
#'
#'  * Default is 0.5.
#' @param resolution Numeric controlling the number of communities found during
#'  clustering where values above one indicate more, smaller communities will
#'  be found.
#'
#'  * Default is 0.1
#' @param mad_coef Numeric scalar specifying the multiplier for the
#'   Median Absolute Deviation (MAD) used to define thresholds for
#'   quality control metrics of damaged clusters, allowing flexibility to make
#'   the detection more or less stringent.
#'
#'    * Default is 3.
#' @param pca_columns String vector containing the names of the metrics used
#'  for principal component analysis.
#'
#'  There are 8 options:
#'   - features
#'   - log.features
#'   - counts
#'   - log.counts
#'   - mt.prop
#'   - rb.prop
#'   - malat1.prop
#'   - malat1
#'
#' Default is c("log.features", "log.counts", "mt.prop", "rb.prop").
#' @param pK Numeric describing how many nearest neighbours are considered
#'  for pANN calculations. Note that pK is expressed as a proportion of the
#'   total cell number.
#'
#'  * Default is 0.3.
#' @param filter_counts Boolean specifying whether the output matrix
#'  should be filtered, returned containing only cells that fall below
#'  the filter threshold. Alternatively, a data frame containing cell
#'  barcodes and their associated label as either 'damaged' or 'cell'
#'  is returned.
#'
#'  * Default is FALSE.
#' @param palette String specifying the three colours that will be used to
#'  create the continuous colour palette for colouring the 'damage_column'.
#'
#'  * Default is a range from purple to red,
#'  c("grey", "#7023FD", "#E60006").
#' @param verbose Boolean specifying whether messages and function progress
#'  should be displayed in the console.
#'
#'  * Default is TRUE.
#' @return Filtered matrix or data frame containing damage labels.
#' @importFrom e1071 skewness
#' @importFrom Seurat CreateSeuratObject NormalizeData ScaleData
#' @importFrom Seurat FindVariableFeatures RunPCA FindNeighbors
#' @importFrom Seurat FindClusters Idents FindMarkers
#' @importFrom stats median
#' @importFrom Matrix colSums Matrix
#' @importFrom RcppHNSW hnsw_knn
#' @importFrom dplyr %>% group_by summarise mutate case_when
#' @importFrom dplyr first pull filter
#' @importFrom ggplot2 ggsave theme element_rect margin
#' @importFrom rlang expr !!! .data
#' @importFrom stats prcomp
#' @importFrom stringr str_extract
#' @export
#' @examples
#' data("test_counts", package = "DamageDetective")
#' test_counts <- test_counts[, 1:100]
#'
#' test <- detect_damage(
#'   count_matrix = test_counts,
#'   generate_plot = FALSE,
#'   seed = 7
#' )
#'
detect_damage <- function(
    count_matrix,
    ribosome_penalty = 1,
    organism = "Hsap",
    annotated_celltypes = FALSE,
    target_damage = c(0.65, 1),
    damage_distribution = "left_skewed",
    distribution_steepness = "steep",
    pca_columns = c("log.features", "log.counts", "mt.prop", "rb.prop"),
    beta_shape_parameters = NULL,
    damage_proportion = 0.15,
    seed = 7,
    pK = 0.3,
    resolution = 0.1,
    mad_coef = 3,
    generate_plot = TRUE,
    display_plot = TRUE,
    palette = c("grey", "#7023FD", "#E60006"),
    filter_threshold = 0.5,
    filter_counts = FALSE,
    verbose = TRUE
) {
  # Data preparation
  .check_inputs(pK, damage_proportion, filter_threshold, count_matrix)
  gene_idx <- get_organism_indices(count_matrix, organism)

  # Cluster cells
  cluster_cells_info <- .cluster_cells(count_matrix, resolution, verbose)
  cluster_cells <- cluster_cells_info$clusters
  seurat_object <- cluster_cells_info$object

  # Check for entirely damaged clusters
  damaged_cluster_cells <- .check_clusters(cluster_cells, seurat_object,
                                           gene_idx, mad_coef, verbose)


  # Generate simulations for defined damage levels
  damaged_matrices <- .simulate_damage(
    count_matrix, cluster_cells, damaged_cluster_cells,
    damage_proportion, ribosome_penalty, organism,
    damage_distribution, target_damage, distribution_steepness,
    seed, verbose
  )
  matrix_combined <- .combine_matrices(damaged_matrices, count_matrix)

  # Compute cell QC metrics
  metadata_stored <- .compute_qc_metrics(matrix_combined, gene_idx)
  metadata_stored_filtered <- metadata_stored[
    !rownames(metadata_stored) %in% damaged_cluster_cells,
  ]

  # Perform PCA and nearest neighbour search
  neighbour_indices <- .perform_pca(
    metadata_stored_filtered, pca_columns, pK, verbose)

  # Compute pANN and scale damage levels
  metadata_plot <- .compile_pANN(
    metadata_stored, neighbour_indices, ranges,
    pK, verbose, damaged_cluster_cells
  )

  # Visualize and filter results
  output <- .finalise_output(
    metadata_plot, count_matrix, filter_counts,
    generate_plot, display_plot, target_damage,
    filter_threshold, palette, damaged_cluster_cells, metadata_stored
  )

  return(output)
}

.check_inputs <- function(pK, damage_proportion,
                          filter_threshold, count_matrix) {
  if (pK > 1 || pK <= 0) {
    stop("pK must lie between 0 and 1.")
  }

  if (damage_proportion > 1) {
    stop("damage_proportion must lie between 0 and 1.")
  }

  if (!inherits(count_matrix, "Matrix")) {
    stop("count_matrix must be a matrix or a sparse matrix.")
  }

  if (!is.numeric(filter_threshold) ||
      filter_threshold <= 0 ||
      filter_threshold > 1
  ) {
    stop("filter_threshold must be a numeric value between 0 and 1.")
  }
}

.cluster_cells <- function(
    count_matrix, resolution, verbose
){
  if (verbose) {
    message("Clustering cells...")
  }

  # Use Seurat's default workflow to identify clusters
  seurat <- suppressWarnings(CreateSeuratObject(counts = count_matrix))
  seurat <- suppressWarnings(NormalizeData(seurat, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE))

  # Retrieve cluster numbers present
  clusters <- table(seurat$seurat_clusters)
  cluster_cells <- list()
  for (cluster in names(clusters)) {
    cells <- subset(seurat, subset = seurat_clusters == cluster)
    cells <- rownames(cells@meta.data)
    cluster_name <- paste0("cluster_", cluster)
    cluster_cells[[cluster_name]] <- cells
  }

  return(list(
    clusters = cluster_cells,
    object = seurat)
  )
}

.check_clusters <- function(
    cluster_cells, seurat_object,
    gene_idx, mad_coef, verbose) {
  if (verbose) message("Checking clusters...")

  Seurat::Idents(seurat_object) <- "seurat_clusters"

  if (length(cluster_cells) <= 1) return(NULL)

  # Compute mitochondrial and ribosomal proportions for each cell
  counts <- Seurat::GetAssayData(seurat_object, slot = "counts")
  mito_genes <- grep(gene_idx$mito_pattern, rownames(seurat_object), value = TRUE)
  ribo_genes <- grep(gene_idx$ribo_pattern, rownames(seurat_object), value = TRUE)

  total_counts <- Matrix::colSums(counts)
  mito_prop <- Matrix::colSums(counts[mito_genes, , drop = FALSE]) / total_counts
  ribo_prop <- Matrix::colSums(counts[ribo_genes, , drop = FALSE]) / total_counts

  cluster_ids <- as.character(Seurat::Idents(seurat_object))

  # Compute medians per cluster
  cluster_medians <- tapply(seq_along(cluster_ids), cluster_ids, function(idx) {
    list(
      mito = median(mito_prop[idx], na.rm = TRUE),
      ribo = median(ribo_prop[idx], na.rm = TRUE)
    )
  })

  mito_medians <- sapply(cluster_medians, function(x) x$mito)
  ribo_medians <- sapply(cluster_medians, function(x) x$ribo)

  # Compute thresholds
  min_mito_median <- min(mito_medians, na.rm = TRUE)
  max_ribo_median <- max(ribo_medians, na.rm = TRUE)
  mito_mad <- stats::mad(mito_prop, constant = 1, na.rm = TRUE)
  ribo_mad <- stats::mad(mito_prop, constant = 1, na.rm = TRUE)
  mito_threshold <- min_mito_median + mad_coef * mito_mad
  ribo_threshold <- max_ribo_median - mad_coef * ribo_mad

  # Identify damaged clusters
  damaged_clusters <- c()
  for (cl in names(cluster_medians)) {
    mito_pass <- mito_medians[cl] >= mito_threshold
    ribo_pass <- ribo_medians[cl] <= ribo_threshold

    if (!is.na(mito_pass) && !is.na(ribo_pass) && mito_pass && ribo_pass) {
      # Match the full cluster name in cluster_cells
      damaged_clusters <- c(damaged_clusters, paste0("cluster_", cl))
    }
  }

  if (length(damaged_clusters) > 0) {
    damaged_cells <- unlist(cluster_cells[damaged_clusters], use.names = FALSE)
  } else {
    damaged_cells <- integer(0)
  }
  return(unique(damaged_cells))
}

.simulate_damage <- function(
    count_matrix, cluster_cells, damaged_cluster_cells,
    damage_proportion, ribosome_penalty, organism,
    damage_distribution, target_damage, distribution_steepness,
    seed, verbose
){
  if (verbose) {
    message("Simulating damage...")
  }

  damaged_matrices <- list()

  # Skip clusters with immediate damaged signature
  for (cluster_name in names(cluster_cells)) {
    cells_to_simulate <- setdiff(cluster_cells[[cluster_name]], damaged_cluster_cells)
    if (length(cells_to_simulate) == 0) {
      next
    }

    matrix_subset <- count_matrix[, cells_to_simulate, drop = FALSE]

    damaged_cells <- simulate_counts(
      count_matrix = matrix_subset,
      damage_proportion = damage_proportion,
      target_damage = target_damage,
      ribosome_penalty = ribosome_penalty,
      damage_distribution = damage_distribution,
      distribution_steepness = distribution_steepness,
      generate_plot = FALSE,
      seed = seed
    )

    barcodes <- damaged_cells$qc_summary %>%
      dplyr::filter(.data$Damaged_Level != 0) %>%
      dplyr::pull(.data$Cell)

    damaged_matrix <- damaged_cells$matrix[, barcodes, drop = FALSE]
    colnames(damaged_matrix) <- paste0(colnames(damaged_matrix), "_", cluster_name)

    damaged_matrices[[cluster_name]] <- damaged_matrix
  }

  return(damaged_matrices)
}


.combine_matrices <- function(damaged_matrices, count_matrix) {
  matrix_updated <- do.call(cbind, damaged_matrices)
  matrix_combined <- cbind(matrix_updated, count_matrix)
  matrix_combined <- matrix_combined[, unique(colnames(matrix_combined))]
  return(matrix_combined)
}

.compute_qc_metrics <- function(matrix_combined, gene_idx, ranges) {

  # Mark damaged cells as 1 and true as 0
  Damage_level <- str_extract(colnames(matrix_combined), "cluster_\\d+")
  Damage_level <- ifelse(is.na(Damage_level), "true", Damage_level)

  # Compute QC metrics
  features <- Matrix::colSums(matrix_combined > 0)
  log.features <- log10(features)
  total_counts <- Matrix::colSums(matrix_combined)
  log.counts <- log10(total_counts)
  mt.prop <- Matrix::colSums(
    matrix_combined[
      grep(gene_idx$mito_pattern, rownames(matrix_combined)), , drop = FALSE
    ]
  ) / total_counts
  epsilon <- 1e-6
  mt.logit <- log((mt.prop + epsilon) / (1 - mt.prop + epsilon))
  rb.prop <- Matrix::colSums(
    matrix_combined[
      grep(gene_idx$ribo_pattern, rownames(matrix_combined)), , drop = FALSE]
  ) / total_counts

  if (gene_idx$MALAT1 %in% seq_len(nrow(matrix_combined))) {
    malat1_expr <- matrix_combined[gene_idx$MALAT1, , drop = FALSE]
    malat1 <- as.vector(malat1_expr)
    malat1.prop <- Matrix::colSums(malat1_expr) / total_counts
  } else {
    # Create zero values if MALAT1 is not in the count matrix
    malat1_expr <- Matrix::Matrix(0, nrow = 1, ncol = ncol(matrix_combined), sparse = TRUE)
    malat1 <- rep(0, ncol(matrix_combined))
    malat1.prop <- rep(0, ncol(matrix_combined))
  }
  malat1.arcsin <- asin(sqrt(malat1.prop))

  # Compile QC dataframe
  metadata_stored <- data.frame(
    features = features,
    log.features = log.features,
    counts = total_counts,
    log.counts = log.counts,
    mt.prop = mt.prop,
    mt.logit = mt.logit,
    rb.prop = rb.prop,
    malat1.prop = malat1.prop,
    malat1.arcsin = malat1.arcsin,
    malat1 = malat1,
    Damage_level = Damage_level,
    row.names = colnames(matrix_combined)
  )
  return(metadata_stored)
}

.perform_pca <- function(metadata_stored, pca_columns, pK, verbose) {

  # Adjust heavily skewed data
  skew_test <- e1071::skewness(metadata_stored$mt.prop) > 1.5
  if (skew_test){
    pca_columns <- c("log.features", "log.counts", "mt.logit", "rb.prop")
  }

  # Test for high expression of MALAT1
  malat_test <- (max(metadata_stored$malat1.prop) > 0.5)
  if (malat_test){
    pca_columns <- c("log.features", "mt.prop", "malat1.prop")
  }

  qc_pcs <- length(pca_columns)
  metadata_pca <- metadata_stored[, pca_columns]
  pca_result <- prcomp(metadata_pca, center = TRUE, scale. = TRUE)
  k <- round(pK * nrow(metadata_stored), 0)
  RcppHNSW::hnsw_knn(pca_result$x[, 1:qc_pcs], k)$idx
}

.compute_pANN <- function(barcode_set, metadata_stored, neighbour_indices, pK) {
  vapply(rownames(metadata_stored), function(cell) {

    # Get the neighbours for the cell (excluding itself)
    index <- which(rownames(metadata_stored) == cell)
    neighbours <- neighbour_indices[index, -1]
    neighbour_barcodes <- rownames(metadata_stored)[neighbours]

    # Compute proportion of neighbours in the barcode set
    k <- round(pK * nrow(metadata_stored), 0)
    sum(neighbour_barcodes %in% barcode_set) / k
  }, FUN.VALUE = numeric(1))
}

.compile_pANN <- function(
    metadata_stored, neighbour_indices, ranges,
    pK, verbose, damaged_cluster_cells
) {
  # Compute proportion of artificial neighbours
  if (verbose) {
    message("Computing pANN...")
  }

  # Isolate cells where damage needs to be found
  metadata_stored <- metadata_stored[!rownames(metadata_stored) %in%
                                       damaged_cluster_cells, ]


  # Isolate columns & ensure cell names are present
  metadata_columns <- c("features", "counts", "mt.prop", "rb.prop",
                        "malat1", "Damage_level")

  metadata_plot <- metadata_stored[, metadata_columns]
  metadata_plot$Cells <- rownames(metadata_stored)
  metadata_plot$Damage_level <- as.character(metadata_plot$Damage_level)

  # Retrieve barcodes based on each artificial cluster
  clusters <- unique(metadata_plot$Damage_level)
  clusters <- clusters[clusters != "true"]

  barcodes <- list()
  for (cluster in clusters) {
    barcodes[[cluster]] <- metadata_plot$Cells[metadata_plot$Damage_level
                                               == cluster]
  }

  # Compute proportion of artificial neighbours for each cluster
  pANN_names <- c()
  for (cluster in clusters) {
    pANN_name <- paste0("pANN_", cluster)
    pANN_names <- append(pANN_name, pANN_names)
    metadata_plot[[pANN_name]] <- .compute_pANN(
      barcodes[[cluster]], metadata_stored, neighbour_indices, pK
    )
  }

  # Normalize according to the population size of each artificial cluster
  n_artificial <- sum(metadata_plot$Damage_level != "true")
  cluster_sizes <- table(metadata_plot$Damage_level)
  cluster_sizes <- subset(cluster_sizes, names(cluster_sizes) != "true")
  cluster_proportions <- cluster_sizes / n_artificial

  # Isolate true cells & define relevant columns
  metadata_plot <- metadata_plot %>%
    dplyr::filter(.data$Damage_level == "true")

  for (cluster in clusters){
    cluster_name <- paste0("pANN_", cluster)
    metadata_plot[[cluster_name]] <- metadata_plot[[cluster_name]] /
      cluster_proportions[[cluster]]
  }

  # Assign the value of the highest pANN
  pann_cols <- grep("^pANN_cluster_", colnames(metadata_plot), value = TRUE)

  if (length(pann_cols) == 1){
    metadata_plot$pANN <- metadata_plot[[pann_cols]]
  } else {
    metadata_plot$pANN <- apply(metadata_plot[, pann_cols], 1, max)
  }

  return(metadata_plot)

}

.finalise_output <- function(metadata_plot, count_matrix, filter_counts,
                             generate_plot, display_plot, target_damage,
                             filter_threshold, palette, damaged_cluster_cells,
                             metadata_stored) {

  # Cells with simulated damage scores
  metadata_plot$DamageDetective <-
    (metadata_plot$pANN - min(metadata_plot$pANN, na.rm = TRUE)) /
    (max(metadata_plot$pANN, na.rm = TRUE) - min(metadata_plot$pANN, na.rm = TRUE))

  simulated_output <- metadata_plot[, c("Cells", "features", "counts", "mt.prop", "rb.prop", "malat1", "DamageDetective")]

  # Damaged cluster cells qc metrics
  if (!is.null(damaged_cluster_cells)) {
    damaged_metadata <- metadata_stored[rownames(metadata_stored) %in% damaged_cluster_cells, ]

    damaged_output <- data.frame(
      Cells = rownames(damaged_metadata),
      features = damaged_metadata$features,
      counts = damaged_metadata$counts,
      mt.prop = damaged_metadata$mt.prop,
      rb.prop = damaged_metadata$rb.prop,
      malat1 = damaged_metadata$malat1,
      DamageDetective = 1
    )

    # Combine & reorder
    metadata_output <- rbind(simulated_output, damaged_output)
    metadata_output <- metadata_output[match(colnames(count_matrix),
                                             metadata_output$Cells), ]
  } else {
    metadata_output <- simulated_output
    metadata_output <- metadata_output[match(colnames(count_matrix),
                                             metadata_output$Cells), ]
  }

  if (filter_counts) {
    filtered_cells <- metadata_output$Cells[metadata_output$DamageDetective <= filter_threshold]
    final_filtered_matrix <- count_matrix[, filtered_cells, drop = FALSE]
    output <- final_filtered_matrix
  } else {
    final_output <- metadata_output[, c("Cells", "DamageDetective")]
    output <- final_output
  }

  # Generate plot if needed
  if (generate_plot) {
    final_plot <- plot_detection_outcome(
      qc_summary = metadata_output,
      target_damage = target_damage,
      palette = palette
    )

    if (display_plot) {
      print(final_plot)
    }

    output <- list(
      output = output,
      plot = final_plot
    )
  }

  return(output)
}
