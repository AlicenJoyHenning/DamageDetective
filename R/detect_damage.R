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
#' @param qc_columns String vector containing the names of the metrics used
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
#' @param parallel_threshold Numeric specifying the maximum number of cells
#'  for which parallel processing will be used. For datasets larger than this
#'  threshold, sequential processing will be used to avoid memory issues with
#'  large data exports to parallel workers.
#'
#'  * Default is 10000.
#' @return Filtered matrix or data frame containing damage labels.
#' @importFrom e1071 skewness
#' @importFrom Seurat CreateSeuratObject NormalizeData ScaleData
#' @importFrom Seurat FindVariableFeatures RunPCA FindNeighbors
#' @importFrom Seurat FindClusters Idents FindMarkers
#' @importFrom Matrix colSums Matrix
#' @importFrom RcppHNSW hnsw_knn
#' @importFrom dplyr %>% group_by summarise mutate case_when
#' @importFrom dplyr first pull filter
#' @importFrom ggplot2 ggsave theme element_rect margin
#' @importFrom rlang expr !!! .data
#' @importFrom stats prcomp
#' @importFrom stringr str_extract
#' @importFrom future.apply future_lapply
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
    target_damage = c(0.7, 1),
    damage_distribution = "left_skewed",
    distribution_steepness = "steep",
    qc_columns = c( "mt.prop"),
    beta_shape_parameters = NULL,
    damage_proportion = 0.15,
    seed = 7,
    resolution = 0.1,
    generate_plot = TRUE,
    display_plot = TRUE,
    palette = c("grey", "#7023FD", "#E60006"),
    filter_threshold = 0.5,
    filter_counts = FALSE,
    parallel_threshold = 3500,
    verbose = TRUE
) {
  # Data preparation
  .check_inputs(filter_threshold, count_matrix)
  gene_idx <- get_organism_indices(count_matrix, organism)

  # Cluster cells
  cluster_cells_info <- .cluster_cells(count_matrix, resolution, verbose)
  cluster_cells <- cluster_cells_info$clusters
  seurat_object <- cluster_cells_info$object

  # Identify entire damaged clusters
  damaged_cluster_cells <- .check_clusters(cluster_cells, seurat_object,
                                           gene_idx, verbose)


  # Generate simulations for defined damage levels
  damaged_matrices <- .simulate_damage(
    count_matrix, cluster_cells, damaged_cluster_cells,
    damage_proportion, ribosome_penalty, organism,
    damage_distribution, target_damage, distribution_steepness, seed, verbose,
    parallel_threshold
  )
  matrix_combined <- .combine_matrices(damaged_matrices, count_matrix)

  # Compute cell QC metrics
  metadata_stored <- .compute_qc_metrics(matrix_combined, gene_idx)
  metadata_stored_filtered <- metadata_stored[
    !rownames(metadata_stored) %in% damaged_cluster_cells,
  ]

  # Perform PCA and nearest neighbour search
  metadata_compared <- .compare_qc(
    metadata_stored = metadata_stored_filtered,
    qc_columns = qc_columns,
    seurat_object = seurat_object,
    verbose = verbose
  )

  # Compute pANN and scale damage levels
  metadata_plot <- .compile_metadata(
    metadata_compared = metadata_compared,
    damaged_cluster_cells = damaged_cluster_cells,
    verbose = verbose
  )

  # Visualize and filter results
  output <- .finalise_output(
    metadata_plot = metadata_plot,
    count_matrix = count_matrix,
    filter_counts = filter_counts,
    generate_plot = generate_plot,
    display_plot = display_plot,
    target_damage = target_damage,
    filter_threshold = filter_threshold,
    palette = palette,
    damaged_cluster_cells = damaged_cluster_cells
  )

  return(output)
}

.check_inputs <- function(filter_threshold, count_matrix) {
  if (!is.numeric(filter_threshold) ||
      filter_threshold <= 0 ||
      filter_threshold > 1
  ) {
    stop("filter_threshold must be a numeric value between 0 and 1.")
  }

  if (!inherits(count_matrix, "Matrix")) {
    stop("count_matrix must be a matrix or a sparse matrix.")
  }
}

.cluster_cells <- function(
    count_matrix, resolution, verbose
){
  if (verbose) {
    message("Clustering cells...")
  }

  # Use Seurat's default workflow to identify clusters
  options(Seurat.threads = parallel::detectCores())
  seurat <- suppressWarnings(CreateSeuratObject(counts = count_matrix))
  seurat <- suppressWarnings(NormalizeData(seurat, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    RunPCA(verbose = FALSE, dims = 1:10) %>%
    FindNeighbors(verbose = FALSE, dims = 1:10) %>%
    FindClusters(resolution = resolution, verbose = FALSE, dims = 1:10))

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

.check_clusters <- function(cluster_cells, seurat_object, gene_idx, verbose) {
  if (verbose) {
    message("Checking clustering signal...")
  }

  Seurat::Idents(seurat_object) <- "seurat_clusters"

  if (length(cluster_cells) == 1) {
    return(NULL)
  } else {

  damaged_cells <- c()

  for (cluster in names(cluster_cells)) {
    cluster_number <- as.numeric(gsub("cluster_", "", cluster))

    markers <- Seurat::FindMarkers(seurat_object,
                           ident.1 = cluster_number,
                           ident.2 = NULL)

    markers <- subset(markers, p_val_adj <= 0.05 & avg_log2FC >= 2)
    markers <- markers[order(-markers$avg_log2FC), ]
    top10_genes <- rownames(head(markers, 10))

    num_matches <- sum(grepl(gene_idx$mito_pattern, top10_genes))
    if (num_matches >= 5) {
      damaged_cells <- c(damaged_cells, cluster_cells[[cluster]])
    }
  }

  return(unique(damaged_cells))

  }
}

.simulate_damage <- function(
    count_matrix, cluster_cells, damaged_cluster_cells,
    damage_proportion, ribosome_penalty, organism,
    damage_distribution, target_damage, distribution_steepness, seed, verbose,
    parallel_threshold
){
  if (verbose) {
    message("Simulating damage...")
  }

  # Filter out clusters with no cells to simulate
  valid_clusters <- names(cluster_cells)[
    vapply(names(cluster_cells), function(cluster_name) {
      cells_to_simulate <- setdiff(cluster_cells[[cluster_name]], damaged_cluster_cells)
      length(cells_to_simulate) > 0
    }, logical(1))
  ]

  if (length(valid_clusters) == 0) {
    return(list())
  }

  # Check if the dataset is large and decide on parallelization strategy
  total_cells <- ncol(count_matrix)
  use_parallel <- total_cells <= parallel_threshold

  if (use_parallel) {
    # For smaller datasets, prepare cluster-specific data to minimize exports
    cluster_data <- lapply(valid_clusters, function(cluster_name) {
      cells_to_simulate <- setdiff(cluster_cells[[cluster_name]], damaged_cluster_cells)
      matrix_subset <- count_matrix[, cells_to_simulate, drop = FALSE]

      # Calculate damage proportion
      cells <- ncol(matrix_subset)
      damage_prop <- if (cells <= 1000) {
        0.7
      } else if (cells <= 2000) {
        0.25
      } else if (cells <= 5000) {
        0.2
      } else {
        0.1
      }

      cluster_seed <- seed + as.numeric(gsub("cluster_", "", cluster_name))

      list(
        matrix_subset = matrix_subset,
        damage_prop = damage_prop,
        cluster_seed = cluster_seed,
        cluster_name = cluster_name
      )
    })
    names(cluster_data) <- valid_clusters

    # Process clusters in parallel using future.apply with pre-prepared data
    damaged_matrices <- future.apply::future_lapply(cluster_data, function(data) {
      damaged_cells <- simulate_counts(
        count_matrix = data$matrix_subset,
        damage_proportion = data$damage_prop,
        target_damage = target_damage,
        ribosome_penalty = ribosome_penalty,
        damage_distribution = damage_distribution,
        distribution_steepness = distribution_steepness,
        generate_plot = FALSE,
        organism = organism,
        seed = data$cluster_seed
      )

      barcodes <- damaged_cells$qc_summary %>%
        dplyr::filter(.data$Damaged_Level != 0) %>%
        dplyr::pull(.data$Cell)

      damaged_matrix <- damaged_cells$matrix[, barcodes, drop = FALSE]
      colnames(damaged_matrix) <- paste0(colnames(damaged_matrix), "_", data$cluster_name)

      return(damaged_matrix)
    }, future.seed = TRUE)
  } else {
    # For large datasets, fall back to sequential processing
    if (verbose) {
      message("Large dataset detected. Using sequential processing to avoid memory issues.")
    }

    damaged_matrices <- list()
    for (cluster_name in valid_clusters) {
      cells_to_simulate <- setdiff(cluster_cells[[cluster_name]], damaged_cluster_cells)
      matrix_subset <- count_matrix[, cells_to_simulate, drop = FALSE]

      # Simulate all as many as feasible without taking too long
      cells <- ncol(matrix_subset)
      damage_prop <- if (cells <= 1000) {
        1
      } else if (cells <= 2000) {
        0.5
      } else if (cells <= 5000) {
        0.2
      } else {
        0.1
      }

      # if (verbose) {
      #   message("Processing cluster ", cluster_name, " with damage proportion: ", damage_prop)
      # }

      # Use unique seed for each cluster to ensure reproducibility
      cluster_seed <- seed + as.numeric(gsub("cluster_", "", cluster_name))

      damaged_cells <- simulate_counts(
        count_matrix = matrix_subset,
        damage_proportion = damage_prop,
        target_damage = target_damage,
        ribosome_penalty = ribosome_penalty,
        damage_distribution = damage_distribution,
        distribution_steepness = distribution_steepness,
        generate_plot = FALSE,
        organism = organism,
        seed = cluster_seed
      )

      barcodes <- damaged_cells$qc_summary %>%
        dplyr::filter(.data$Damaged_Level != 0) %>%
        dplyr::pull(.data$Cell)

      damaged_matrix <- damaged_cells$matrix[, barcodes, drop = FALSE]
      colnames(damaged_matrix) <- paste0(colnames(damaged_matrix), "_", cluster_name)

      damaged_matrices[[cluster_name]] <- damaged_matrix
    }
  }

  # Name the list elements (if not already named)
  if (use_parallel) {
    names(damaged_matrices) <- valid_clusters
  }

  return(damaged_matrices)
}

.combine_matrices <- function(damaged_matrices, count_matrix) {
  matrix_updated <- do.call(cbind, damaged_matrices)
  matrix_combined <- cbind(matrix_updated, count_matrix)
  #matrix_combined <- matrix_combined[, unique(colnames(matrix_combined))]
  return(matrix_combined)
}

.compute_qc_metrics <- function(matrix_combined, gene_idx) {

  # Mark damaged cells as 1 and true as 0
  Damage_level <- str_extract(colnames(matrix_combined), "cluster_\\d+")
  Damage_level <- ifelse(is.na(Damage_level), "true", Damage_level) # true cells identified

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

  # Compile QC data frame
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

.compare_qc <- function(metadata_stored, qc_columns, seurat_object, verbose) {
  if (verbose){
    message("Comparing QC metrics...")
  }

  # Test for high expression of MALAT1
  malat_test <- (max(metadata_stored$malat1.prop) > 0.5)
  if (malat_test){
    qc_columns <- append("malat1.prop", qc_columns)
  }

  # Add original cluster label to true cells
  metadata_stored$true_cluster <- suppressWarnings(ifelse(
    row.names(metadata_stored) %in% colnames(seurat_object),
    seurat_object$seurat_clusters[match(row.names(metadata_stored),
                                        colnames(seurat_object))],
    "artificial"
  ))
  metadata_stored$true_cluster <- suppressWarnings(paste0("cluster_",
                 as.numeric(metadata_stored$true_cluster) - 1))

  # Prepare empty column for storage
  metadata_stored$QC_position <- 0

  # Compare distribution in each cluster for matching true & artificial cells
  clusters <- unique(metadata_stored$Damage_level)
  clusters <- setdiff(clusters, "true")

  for (cluster in clusters){
    QC_position <- list()

    # Find means of QC metric of interest for artificial & true populations
    for(qc_column in qc_columns){
      # Find the means for each group
      true_QC_mean <- mean(
        metadata_stored[
          metadata_stored$true_cluster == cluster,
          qc_column, drop = TRUE]
      )

      artificial_QC_mean <- mean(
        metadata_stored[
          metadata_stored$Damage_level == cluster,
          qc_column, drop = TRUE]
      )

      # Compute the position of a true cell value relative to the two means
      true_values <- metadata_stored[
        metadata_stored$true_cluster == cluster,
        qc_column, drop = TRUE]

      # Direct comparison adjusted to meet assumptions of damage
      rescaled_values <- numeric(length(true_values))
      if (qc_column %in% c("mt.prop", "malat1.prop")){
        rescaled_values[true_values <= true_QC_mean] <- 0
        rescaled_values[true_values >= artificial_QC_mean] <- 1
        in_range <- which(true_values > true_QC_mean & true_values
                          < artificial_QC_mean)
        rescaled_values[in_range] <- (true_values[in_range] - true_QC_mean) /
          (artificial_QC_mean - true_QC_mean)
      } else {
        rescaled_values[true_values >= true_QC_mean] <- 0
        rescaled_values[true_values <= artificial_QC_mean] <- 1
        in_range <- which(true_values > artificial_QC_mean & true_values < true_QC_mean)
        rescaled_values[in_range] <- (true_values[in_range] - true_QC_mean) /
          (artificial_QC_mean - true_QC_mean)
      }
      # Store outcome for each metric
      QC_position[[qc_column]] <- rescaled_values
    }

    # Find mean across QC metrics
    QC_position_df <- as.data.frame(QC_position)
    QC_position_values <- rowMeans(QC_position_df, na.rm = TRUE)
    metadata_stored[metadata_stored$true_cluster == cluster, ]$QC_position <-
      QC_position_values
  }
  return(metadata_stored)
}


.compile_metadata <- function(
    metadata_compared, damaged_cluster_cells, verbose
) {
  # Compute proportion of artificial neighbours
  if (verbose) {
    message("Compiling output...")
  }

  # If present, mark damaged cluster highly
  metadata_compared$QC_position <- ifelse(
    metadata_compared$Damage_level %in% damaged_cluster_cells, 1,
    metadata_compared$QC_position)

  # Isolate true cells
  metadata_compared <- subset(metadata_compared,
                            metadata_compared$Damage_level == "true")
  metadata_compared$Damage_level <- NULL

  return(metadata_compared)

}

.finalise_output <- function(metadata_plot, count_matrix, filter_counts,
                             generate_plot, display_plot, target_damage,
                             filter_threshold, palette, damaged_cluster_cells) {

  # Clean the output columns
  metadata_plot$DamageDetective <- metadata_plot$QC_position
  metadata_plot <- metadata_plot[,
    c("features", "counts", "mt.prop", "rb.prop", "malat1", "DamageDetective"),
    drop = FALSE]


  if (filter_counts) {
    filtered_data <- subset(
      metadata_plot,
      metadata_plot$DamageDetective <= filter_threshold)
    filtered_cells <- rownames(filtered_data)
    final_filtered_matrix <- count_matrix[, filtered_cells, drop = FALSE]
    output <- final_filtered_matrix

  } else {
    final_output <- metadata_plot[, c("DamageDetective"), drop = FALSE]
    row.names(final_output) <- row.names(metadata_plot)
    output <- final_output
  }

  # Generate plot if specified
  target_damage <- 0.5
  if (generate_plot) {
    final_plot <- plot_detection_outcome(
      qc_summary = metadata_plot,
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

