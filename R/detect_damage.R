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
#' Before simulation, cells are clustered on their mitochondrial and ribosomal
#' gene expression and each cluster is labelled as 'viable', 'divergent' or
#' 'depleted' by comparing its median QC metrics to those of the cluster with
#' the lowest non-trivial mitochondrial percentage. Only viable cells are used
#' for simulation. Cells in depleted or divergent clusters are excluded from
#' simulation and the pANN search, and are assigned the maximum damage score;
#' cells in divergent clusters are additionally labelled 'divergent' in the
#' output.
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
#' @param annotated_celltypes Data frame containing cell type annotations,
#' with cell barcodes as row names and the cell type of each cell stored in a
#' column named 'annotation'. If provided, simulation is performed within each
#' cell type rather than within each cluster identified by the cluster check.
#'
#' * Default is NULL.
#' @param filter_threshold Numeric specifying the proportion of RNA loss
#' above which a cell should be considered damaged.
#'
#' * Default is 0.5.
#' @param resolution Numeric controlling the number of communities found during
#' clustering where values above one indicate more, smaller communities will
#' be found.
#'
#' * Default is 0.1
#' @param mad_coef Numeric. No longer used by the cluster check (retained so
#' existing calls do not break).
#'
#' * Default is 3.
#' @param pca_columns String vector containing the names of the metrics used
#' for principal component analysis.
#'
#' There are 8 options:
#' - features
#' - log.features
#' - counts
#' - log.counts
#' - mt.prop
#' - rb.prop
#' - malat1.prop
#' - malat1
#'
#' Default is c("log.features", "log.counts", "mt.prop", "rb.prop").
#' @param pK Numeric describing how many nearest neighbours are considered
#' for pANN calculations. Note that pK is expressed as a proportion of the
#' total cell number.
#'
#' * Default is 0.3.
#' @param filter_counts Boolean specifying whether the output matrix
#' should be filtered, returned containing only cells that fall below
#' the filter threshold. Alternatively, a data frame containing cell
#' barcodes and their associated label as either 'damaged' or 'cell'
#' is returned.
#'
#' * Default is FALSE.
#' @param palette String specifying the three colours that will be used to
#' create the continuous colour palette for colouring the 'damage_column'.
#'
#' * Default is a range from purple to red,
#' c("grey", "#7023FD", "#E60006").
#' @param verbose Boolean specifying whether messages and function progress
#' should be displayed in the console.
#'
#' * Default is TRUE.
#' @return Filtered matrix or data frame containing damage labels ('cell',
#' 'damaged' or 'divergent'). If
#' `generate_plot = TRUE`, a list containing `output`, `plot` (the detection
#' outcome), `pann_plot` (the pANN threshold density, only when the threshold
#' is determined automatically) and `presim_plot` (the pre-simulation cluster
#' check).
#' @importFrom e1071 skewness
#' @importFrom Seurat CreateSeuratObject NormalizeData ScaleData
#' @importFrom Seurat FindVariableFeatures RunPCA FindNeighbors
#' @importFrom Seurat VariableFeatures FindClusters Idents FindMarkers
#' @importFrom Seurat GetAssayData FeaturePlot DimPlot
#' @importFrom stats median setNames
#' @importFrom Matrix colSums Matrix
#' @importFrom RcppHNSW hnsw_knn
#' @importFrom dplyr %>% group_by summarise mutate case_when
#' @importFrom dplyr first pull filter slice_min
#' @importFrom ggplot2 ggsave theme element_rect element_blank margin
#' @importFrom patchwork wrap_plots
#' @importFrom rlang expr !!! .data
#' @importFrom stats density prcomp
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
    ribosome_penalty = 0.5,
    organism = "Hsap",
    annotated_celltypes = NULL,
    target_damage = c(0.65, 1),
    damage_distribution = "left_skewed",
    distribution_steepness = "steep",
    pca_columns = c("log.features", "log.counts", "mt.prop", "rb.prop"),
    beta_shape_parameters = NULL,
    damage_proportion = 0.2,
    seed = 7,
    pK = 0.25,
    resolution = 0.2,
    mad_coef = 3,
    generate_plot = TRUE,
    display_plot = TRUE,
    palette = c("grey", "#7023FD", "#E60006"),
    filter_threshold = NA,
    filter_counts = FALSE,
    verbose = TRUE
) {

  # Data preparation
  if (isFALSE(annotated_celltypes)) annotated_celltypes <- NULL
  .check_inputs(pK, damage_proportion, filter_threshold, count_matrix,
                annotated_celltypes)
  gene_idx <- get_organism_indices(count_matrix, organism)

  # Cluster cells and check clusters for entirely damaged populations
  cluster_info <- .check_clusters(
    count_matrix, resolution, gene_idx, generate_plot, verbose
  )
  cluster_cells         <- cluster_info$clusters
  cells_to_simulate     <- cluster_info$cells_to_simulate
  damaged_cluster_cells <- cluster_info$damaged_cells
  divergent_cells       <- cluster_info$divergent_cells
  presim_plot           <- cluster_info$presim_plot

  # Generate simulations for defined damage levels
  damaged_matrices <- .simulate_damage(
    count_matrix, cluster_cells, cells_to_simulate, annotated_celltypes,
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
    filter_threshold, palette, damaged_cluster_cells, metadata_stored,
    presim_plot, divergent_cells, verbose
  )

  return(output)
}


.check_inputs <- function(pK, damage_proportion,
                          filter_threshold, count_matrix,
                          annotated_celltypes = NULL) {
  if (pK > 1 || pK <= 0) {
    stop("pK must lie between 0 and 1.")
  }
  if (damage_proportion > 1) {
    stop("damage_proportion must lie between 0 and 1.")
  }
  if (!is.na(filter_threshold) &&
      (!is.numeric(filter_threshold) ||
       filter_threshold <= 0 ||
       filter_threshold > 1)
  ) {
    stop("filter_threshold must be a numeric value between 0 and 1.")
  }
  if (!is.null(annotated_celltypes)) {
    if (!is.data.frame(annotated_celltypes) ||
        !"annotation" %in% colnames(annotated_celltypes)) {
      stop("annotated_celltypes must be a data frame with cell barcodes as ",
           "row names and an 'annotation' column.")
    }
    if (!all(colnames(count_matrix) %in% rownames(annotated_celltypes))) {
      warning("Some cells in count_matrix have no annotation and will not ",
              "be used for simulation.")
    }
  }
}


.check_clusters <- function(
    count_matrix, resolution, gene_idx, generate_plot, verbose
){
  if (verbose) {
    message("Clustering cells...")
  }

  # Create Seurat object
  seurat <- suppressWarnings(Seurat::CreateSeuratObject(counts = count_matrix))

  # Isolate genes of interest
  mito_genes <- grep(gene_idx$mito_pattern, rownames(seurat), value = TRUE)
  ribo_genes <- grep(gene_idx$ribo_pattern, rownames(seurat), value = TRUE)
  mt_rb_genes <- unique(c(mito_genes, ribo_genes))

  # Compute QC stats per cell
  counts <- Seurat::GetAssayData(seurat, layer = "counts")
  total_counts <- Matrix::colSums(counts)
  mt_exp <- Matrix::colSums(counts[mito_genes, , drop = FALSE])
  pctMT <- (mt_exp / total_counts) * 100
  pctMT[is.na(pctMT)] <- 0
  seurat$pctMT <- as.numeric(pctMT)
  seurat$mt_exp <- as.numeric(mt_exp)

  # Reduce dimensions to mt-rb
  seurat_reduced <- subset(seurat, features = mt_rb_genes)

  # Standard workflow with reduced memory footprint
  seurat_reduced <- suppressWarnings(
    Seurat::NormalizeData(seurat_reduced, verbose = FALSE))
  seurat_reduced <- suppressWarnings(
    Seurat::FindVariableFeatures(seurat_reduced, verbose = FALSE))
  seurat_reduced <- suppressWarnings(
    Seurat::ScaleData(seurat_reduced, verbose = FALSE))
  seurat_reduced <- suppressWarnings(
    Seurat::RunPCA(seurat_reduced, verbose = FALSE))
  seurat_reduced <- suppressWarnings(
    Seurat::FindNeighbors(seurat_reduced, dims = 1:20, verbose = FALSE))
  seurat_reduced <- suppressWarnings(
    Seurat::FindClusters(seurat_reduced,
                         resolution = resolution,
                         verbose = FALSE))
  seurat_reduced <- suppressWarnings(
    Seurat::RunUMAP(seurat_reduced, dims = 1:20, verbose = FALSE))
  

  seurat_reduced$nFeature_RNA <- seurat$nFeature_RNA[
    colnames(seurat_reduced)
  ]

  if (verbose) message("Checking clusters...")

  # Compute cluster QC
  cluster_QC <- seurat_reduced@meta.data %>%
    dplyr::group_by(.data$seurat_clusters) %>%
    dplyr::summarise(
      pct_MT   = stats::median(.data$pctMT, na.rm = TRUE),
      mt_exp   = stats::median(.data$mt_exp, na.rm = TRUE),
      features = stats::median(.data$nFeature_RNA, na.rm = TRUE)
    )

  # Isolate cluster with lowest non-trivial pctMT
  highQC_cluster <- cluster_QC %>%
    dplyr::filter(.data$pct_MT >= 1) %>%
    dplyr::slice_min(.data$pct_MT, n = 1, with_ties = FALSE)

  if (nrow(highQC_cluster) == 0) {
    # No cluster has a non-trivial pctMT, so there is nothing to compare to
    cluster_QC$QC_label <- "viable"
  } else {

    # Define cluster-level thresholds
    pctMT_threshold  <- highQC_cluster$pct_MT * 2
    mt_exp_threshold <- highQC_cluster$mt_exp * 2
    feature_threshold <- highQC_cluster$features * 0.5

    # Classify clusters
    cluster_QC <- cluster_QC %>%
      dplyr::mutate(
        QC_label = dplyr::case_when(
          # Passive depletion
          pct_MT >= 50 & ((features < 100) |
                            (features <= feature_threshold)) ~ "depleted",
          # Active depletion
          pct_MT >= pctMT_threshold &
            ((features <= feature_threshold) |
               (abs(features - highQC_cluster$features) >
                  abs(features - feature_threshold))) ~ "depleted",
          # Active dysregulation
          (pct_MT >= pctMT_threshold) |
            ((abs(pct_MT - highQC_cluster$pct_MT) >
                abs(pct_MT - pctMT_threshold)) &
               mt_exp >= mt_exp_threshold &
               features >= feature_threshold) ~ "divergent",
          TRUE ~ "viable"
        )
      )
  }

  # Assign cluster QC_labels
  cluster_to_label <- stats::setNames(
    cluster_QC$QC_label,
    as.character(cluster_QC$seurat_clusters)
  )
  seurat_reduced$QC_label <- unname(cluster_to_label[as.character(
    seurat_reduced$seurat_clusters)]
  )

  # Visualise
  presim_plot <- NULL
  if (generate_plot){
    plot_theme <- ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black")
    )
    pctMT_view <- Seurat::FeaturePlot(
      seurat_reduced, features = "pctMT") + plot_theme
    feature_view <- Seurat::FeaturePlot(
      seurat_reduced, features = "nFeature_RNA") + plot_theme
    labelled_view <- Seurat::DimPlot(
      seurat_reduced, group.by = "QC_label") + plot_theme
    presim_plot <- patchwork::wrap_plots(
      pctMT_view, feature_view, labelled_view, nrow = 1)
  }

  # Transfer labels to seurat
  seurat$cluster_QC <- seurat_reduced$QC_label[match(
    colnames(seurat), colnames(seurat_reduced)
  )]
  seurat$seurat_clusters <- seurat_reduced$seurat_clusters[match(
    colnames(seurat), colnames(seurat_reduced)
  )]

  # Isolate cells to simulate damage in (viable) and cells to exclude from
  # simulation (depleted & divergent), keeping track of the divergent cells
  cells_to_simulate <- colnames(seurat)[seurat$cluster_QC == "viable"]
  damaged_cells <- colnames(seurat)[seurat$cluster_QC != "viable"]
  divergent_cells <- colnames(seurat)[seurat$cluster_QC == "divergent"]

  # Retrieve cluster membership efficiently
  clusters <- split(rownames(seurat@meta.data), seurat$seurat_clusters)
  names(clusters) <- paste0("cluster_", names(clusters))

  return(list(
    presim_plot = presim_plot,
    object = seurat,
    clusters = clusters,
    cells_to_simulate = unique(cells_to_simulate),
    damaged_cells = unique(damaged_cells),
    divergent_cells = unique(divergent_cells)
  ))
}

.simulate_damage <- function(
    count_matrix, cluster_cells, cells_to_simulate, annotated_celltypes,
    damage_proportion, ribosome_penalty, organism,
    damage_distribution, target_damage, distribution_steepness,
    seed, verbose
){
  if (verbose) {
    message("Simulating damage...")
  }

  # Redefine groups if annotations provided
  if (!is.null(annotated_celltypes)) {
    cluster_cells <- split(
      colnames(count_matrix),
      annotated_celltypes[colnames(count_matrix), "annotation"]
    )
    names(cluster_cells) <- paste0("cell_type_", names(cluster_cells))
  }

  damaged_matrices <- list()

  # Only simulate in viable cells
  for (cluster_name in names(cluster_cells)) {

    cells_in_cluster <- intersect(cluster_cells[[cluster_name]], cells_to_simulate)

    if (length(cells_in_cluster) == 0) {
      next
    }

    matrix_subset <- count_matrix[, cells_in_cluster, drop = FALSE]

    damaged_cells <- simulate_counts(
      count_matrix = matrix_subset,
      damage_proportion = damage_proportion,
      target_damage = target_damage,
      ribosome_penalty = ribosome_penalty,
      organism = organism,
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
  Damage_level <- str_extract(
    colnames(matrix_combined), "(cluster_\\d+|cell_type_.+)$")
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
  skew_val <- e1071::skewness(metadata_stored$mt.prop, na.rm = TRUE)
  skew_test <- !is.na(skew_val) && skew_val > 1.5

  if (skew_test){
    pca_columns <- c("log.features", "log.counts", "mt.logit", "rb.prop")
  }

  # Test for high expression of MALAT1
  malat_test <- (max(metadata_stored$malat1.prop) >= 0.5)

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
  pann_cols <- paste0("pANN_", clusters)

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
                             metadata_stored, presim_plot, divergent_cells,
                             verbose) {

  # Cells with simulated damage scores
  metadata_plot$DamageDetective <-
    (metadata_plot$pANN - min(metadata_plot$pANN, na.rm = TRUE)) /
    (max(metadata_plot$pANN, na.rm = TRUE) - min(metadata_plot$pANN, na.rm = TRUE))

  simulated_output <- metadata_plot[, c("Cells", "features", "counts", "mt.prop", "rb.prop", "malat1", "DamageDetective")]

  # Damaged cluster cells qc metrics
  if (!is.null(damaged_cluster_cells)) {

    damaged_metadata <- metadata_stored[rownames(metadata_stored) %in% damaged_cluster_cells, ]

    if (nrow(damaged_metadata) == 0) {
      metadata_output <- simulated_output
      metadata_output <- metadata_output[match(colnames(count_matrix),
                                               metadata_output$Cells), ]
    } else {
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
    }
  } else {
    metadata_output <- simulated_output
    metadata_output <- metadata_output[match(colnames(count_matrix),
                                             metadata_output$Cells), ]
  }

  # Dynamically find pANN threshold if not specified
  pann_plot <- NULL

  if (is.na(filter_threshold)){

    # Compute pANN density
    dens <- density(metadata_output$DamageDetective, from = 0, to = 1, n = 1000)
    y <- dens$y
    x <- dens$x
    dy <- diff(y)

    # Identify peaks
    peaks_idx <- which(diff(sign(dy)) == -2) + 1

    # Handle case with no peaks
    if (length(peaks_idx) < 1) {
      warning("No pANN peaks detected.")
      filter_threshold <- 0.5
      valley_idx <- which.min(y)
    } else {
      # Find local minimum after first peak
      first_peak_idx <- peaks_idx[1]
      search_end <- if (length(peaks_idx) > 1) peaks_idx[2] else length(y)
      valley_idx <- which.min(y[first_peak_idx:search_end]) + first_peak_idx - 1
      filter_threshold <- x[valley_idx]
    }

    if (generate_plot){
      plot_data <- data.frame(x = x, y = y)

      pann_plot <- ggplot2::ggplot(plot_data,
                                   ggplot2::aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_vline(xintercept = filter_threshold,
                            color = "red", linetype = "dashed",
                            linewidth = 0.8) +
        ggplot2::theme_classic(base_size = 17) +
        ggplot2::theme(
          axis.ticks = ggplot2::element_blank(),
          axis.title = ggplot2::element_text(size = 15, face = "bold"),
          plot.title = ggplot2::element_text(size = 17, face = "bold", hjust = 0.5),
          panel.background = ggplot2::element_rect(fill = "#F1F1F1", colour = "black")
        ) +
        ggplot2::labs(
          title = "pANN threshold determination",
          x = "pANN",
          y = "Density"
        )
    }

    if (verbose){
      message(" ... ", round(filter_threshold, 2), " suggested threshold.")
    }
  }

  # Perform filtering if specified
  if (filter_counts) {

    # Perform filtering
    filtered_cells <- metadata_output$Cells[
      metadata_output$DamageDetective <= filter_threshold &
        !metadata_output$Cells %in% divergent_cells]
    final_filtered_matrix <- count_matrix[, filtered_cells, drop = FALSE]
    output <- final_filtered_matrix

  } else {
    final_output <- metadata_output[, c("Cells", "DamageDetective")]
    colnames(final_output)[2] <- "pANN"   # numeric score

    final_output$DamageDetective <- ifelse(
      final_output$pANN <= filter_threshold, "cell", "damaged"
    )

    # Cells from divergent clusters are labelled separately
    final_output$DamageDetective[final_output$Cells %in% divergent_cells] <- "divergent"

    final_output <- final_output[, c("Cells", "pANN", "DamageDetective")]
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
      if (!is.null(presim_plot)) print(presim_plot)
      if (!is.null(pann_plot)) print(pann_plot)
      print(final_plot)
    }

    output <- list(
      output = output,
      plot = final_plot,
      pann_plot = pann_plot,
      presim_plot = presim_plot
    )
  }

  if (verbose){
    message("Completed.")
  }

  return(output)
}