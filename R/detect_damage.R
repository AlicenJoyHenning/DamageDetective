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
#'  * Default is 0.75.
#' @param damage_levels Numeric specifying the number of distinct sets of
#'  artificial damaged cells simulated, each with a defined range of loss.
#'  Default ptions include,
#'
#'  * 3 : c(0.00001, 0.08), c(0.1, 0.4), c(0.5, 0.9)
#'  * 5 : c(0.00001, 0.08), c(0.1, 0.3), c(0.3, 0.5), c(0.5, 0.7), c(0.7, 0.9)
#'  * 7 : c(0.00001, 0.08), c(0.1, 0.3), c(0.3, 0.4), c(0.4, 0.5), c(0.5, 0.7),
#'  c(0.7, 0.9), c(0.9, 0.99999).
#'
#'  A user can also provide a list specifying sets with their own ranges of
#'  loss,
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
#' @importFrom Matrix colSums Matrix
#' @importFrom RcppHNSW hnsw_knn
#' @importFrom dplyr %>% group_by summarise mutate case_when first pull
#' @importFrom ggplot2 ggsave theme element_rect margin
#' @importFrom rlang expr !!! .data
#' @importFrom stats prcomp
#' @export
#' @examples
#' data("test_counts", package = "DamageDetective")
#'
#' test <- detect_damage(
#'   count_matrix = test_counts,
#'   ribosome_penalty = 0.001,
#'   damage_levels = 3,
#'   damage_proportion = 0.1,
#'   generate_plot = FALSE,
#'   seed = 7
#' )
detect_damage <- function(
    count_matrix,
    ribosome_penalty = 0.01,
    organism = "Hsap",
    annotated_celltypes = FALSE,
    target_damage = c(0.1, 0.8),
    damage_distribution = "right_skewed",
    distribution_steepness = "moderate",
    pca_columns = c("log.features", "log.counts", "mt.prop", "rb.prop"),
    beta_shape_parameters = NULL,
    damage_levels = 5,
    damage_proportion = 0.15,
    seed = 7,
    mito_quantile = 0.75,
    kN = NULL,
    generate_plot = TRUE,
    display_plot = TRUE,
    palette = c("grey", "#7023FD", "#E60006"),
    filter_threshold = 0.7,
    filter_counts = FALSE,
    verbose = TRUE
) {
  # Data preparation
  .check_inputs(filter_threshold, mito_quantile,
                damage_levels, count_matrix, kN)
  gene_idx <- get_organism_indices(count_matrix, organism)
  mito_ribo_data <- .compute_mito_ribo(count_matrix, gene_idx)
  filtered_cells <- .filter_cells_by_mito(mito_ribo_data, mito_quantile)
  matrix_filtered <- count_matrix[, filtered_cells, drop = FALSE]

  # Generate simulations for defined damage levels
  ranges <- .get_damage_ranges(damage_levels)
  damaged_matrices <- .simulate_damage(
    matrix_filtered, ranges, damage_proportion, ribosome_penalty,
    damage_distribution, distribution_steepness, seed, verbose
  )
  matrix_combined <- .combine_matrices(damaged_matrices, count_matrix)

  # Compute quality control metrics
  metadata_stored <- .compute_qc_metrics(matrix_combined, gene_idx, ranges)

  # Perform PCA and nearest neighbour search
  kN <- .find_knn(count_matrix, kN)
  neighbour_indices <- .perform_pca(metadata_stored, pca_columns, kN)

  # Compute pANN and scale damage levels
  metadata_plot <- .compile_pANN(
    metadata_stored, neighbour_indices, ranges, kN, verbose
  )
  metadata_plot <- .scale_damage(metadata_plot, ranges, filter_threshold)

  # Visualize and filter results
  output <- .finalise_output(metadata_plot, count_matrix, filter_counts,
                             generate_plot, display_plot, target_damage,
                             filter_threshold, palette)

  return(output)
}

.check_inputs <- function(filter_threshold, mito_quantile, damage_levels,
                          count_matrix, kN
) {
  if (!is.numeric(filter_threshold) ||
      filter_threshold <= 0 ||
      filter_threshold > 1
  ) {
    stop("filter_threshold must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(mito_quantile) || mito_quantile <= 0 || mito_quantile > 1) {
    stop("mito_quantile must be a numeric value
         0 and 1.")
  }
  if (!is.numeric(damage_levels) && !is.list(damage_levels)) {
    stop("damage_levels must be a numeric value or a list of ranges.")
  }
  if (!is.matrix(count_matrix) && !inherits(count_matrix, "Matrix")) {
    stop("count_matrix must be a matrix or a sparse matrix.")
  }
  if (!is.null(kN) && (!is.numeric(kN) || kN >= 0 ||
                       kN <= dim(count_matrix)[2])) {
    stop("kN must be a positive numeric value less than the number of
         cells present.")
  }
}

.compute_mito_ribo <- function(count_matrix, gene_idx) {
  mito <- colSums(count_matrix[gene_idx$mito_idx, , drop = FALSE]) /
    colSums(count_matrix)
  ribo <- colSums(count_matrix[gene_idx$ribo_idx, , drop = FALSE]) /
    colSums(count_matrix)
  mito_ribo_data <- cbind(mito, ribo)
  colnames(mito_ribo_data) <- c("mito", "ribo")
  rownames(mito_ribo_data) <- colnames(count_matrix)

  if (mean(mito) == 0){
    stop("Please ensure `organism` is correct.")
  }

  return(mito_ribo_data)
}

.filter_cells_by_mito <- function(mito_ribo_data, mito_quantile) {
  mito_top <- quantile(
    mito_ribo_data[, "mito"],
    probs = mito_quantile,
    na.rm = TRUE
  )
  filtered_cells <- rownames(mito_ribo_data)[mito_ribo_data[, "mito"] <=
                                               mito_top]
  return(filtered_cells)
}

.get_damage_ranges <- function(damage_levels) {
  if (is.list(damage_levels)) {
    return(damage_levels)
  }
  if (damage_levels == 3) {
    return(list(
      pANN_10 = c(0.00001, 0.08),
      pANN_40 = c(0.1, 0.4),
      pANN_90 = c(0.5, 0.9)
    ))
  }
  if (damage_levels == 5) {
    return(list(
      pANN_10 = c(0.00001, 0.08),
      pANN_30 = c(0.1, 0.3),
      pANN_50 = c(0.3, 0.5),
      pANN_70 = c(0.5, 0.7),
      pANN_90 = c(0.7, 0.9)
    ))
  }
  if (damage_levels == 7) {
    return(list(
      pANN_10 = c(0.00001, 0.08),
      pANN_30 = c(0.1, 0.3),
      pANN_40 = c(0.3, 0.4),
      pANN_50 = c(0.4, 0.5),
      pANN_70 = c(0.5, 0.7),
      pANN_90 = c(0.7, 0.9),
      pANN_100 = c(0.9, 0.99999)
    ))
  }
}

.simulate_damage <- function(
    matrix_filtered, ranges, damage_proportion, ribosome_penalty,
    damage_distribution, distribution_steepness, seed, verbose
){
  damaged_matrices <- list()
  damage_ranges <- c()

  for (i in 1:length(ranges)) {
    # Initialize list to store damaged matrices
    damaged_matrices <- list()
    damage_ranges <- c()

    # Iterate through ranges and simulate damage
    for (i in 1:length(ranges)) {
      level <- names(ranges)[i]
      damage_range <- ranges[[i]]
      damage_ranges <- append(damage_range, damage_ranges)

      if (verbose) {
        message("Simulating ",
                damage_range[[1]], " and ", damage_range[[2]], " RNA loss...")
      }

      # Simulate counts
      damaged_cells <- simulate_counts(
        count_matrix = matrix_filtered,
        damage_proportion = damage_proportion,
        target_damage = damage_range,
        ribosome_penalty = ribosome_penalty,
        damage_distribution = damage_distribution,
        distribution_steepness = distribution_steepness,
        generate_plot = FALSE,
        seed = seed
      )

      # Extract barcodes of damaged cells
      barcodes <- damaged_cells$qc_summary %>%
        dplyr::filter(.data$Damaged_Level != 0) %>%
        dplyr::pull(.data$Cell)

      damaged_matrix <- damaged_cells$matrix[, barcodes]

      # Append the level of damage to the barcodes
      colnames(damaged_matrix) <- paste0(
        colnames(damaged_matrix), "_", gsub("pANN_", "", level)
      )

      # Store the damaged matrix in the list
      damaged_matrices[[level]] <- damaged_matrix

    }
    return(damaged_matrices)
  }
}

.combine_matrices <- function(damaged_matrices, count_matrix) {
  matrix_updated <- do.call(cbind, damaged_matrices)
  matrix_combined <- cbind(matrix_updated, count_matrix)
  matrix_combined <- matrix_combined[, unique(colnames(matrix_combined))]
  return(matrix_combined)
}

.compute_qc_metrics <- function(matrix_combined, gene_idx, ranges) {

  # Extract Damage_level from cell names
  damage_numbers <- sub(".*_", "", names(ranges))
  cell_names <- colnames(matrix_combined)
  Damage_level <- sub(".*_", "", cell_names)
  Damage_level <- ifelse(!Damage_level %in% damage_numbers, 0, Damage_level)

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
  rb.prop <- Matrix::colSums(
    matrix_combined[
      grep(gene_idx$ribo_pattern, rownames(matrix_combined)), , drop = FALSE]
  ) / total_counts
  malat1_expr <- matrix_combined[gene_idx$MALAT1, , drop = FALSE]
  malat1 <- as.vector(malat1_expr)
  malat1.prop <- Matrix::colSums(
    matrix_combined[
      gene_idx$MALAT1, , drop = FALSE]
  ) / total_counts

  # Compile QC dataframe
  metadata_stored <- data.frame(
    features = features,
    log.features = log.features,
    counts = total_counts,
    log.counts = log.counts,
    mt.prop = mt.prop,
    rb.prop = rb.prop,
    malat1.prop = malat1.prop,
    malat1 = malat1,
    Damage_level = Damage_level,
    row.names = colnames(matrix_combined)
  )
  return(metadata_stored)
}

.find_knn <- function(count_matrix, kN) {
  if (is.null(kN)) {
    kN <- round(dim(count_matrix)[2] / 5, 0)
  }
  return(kN)
}

.perform_pca <- function(metadata_stored, pca_columns, kN) {
  # pca_columns <- c("log.features", "log.counts", "mt.prop",
  # "rb.prop", "malat1")
  qc_pcs <- length(pca_columns)
  metadata_pca <- metadata_stored[, pca_columns]
  pca_result <- prcomp(metadata_pca, center = TRUE, scale. = TRUE)
  RcppHNSW::hnsw_knn(pca_result$x[, 1:qc_pcs], k = kN)$idx
}

.compute_pANN <- function(barcode_set, metadata_stored, neighbour_indices, kN) {
  vapply(rownames(metadata_stored), function(cell) {

    # Get the neighbours for the cell (excluding itself)
    index <- which(rownames(metadata_stored) == cell)
    neighbours <- neighbour_indices[index, -1]
    neighbour_barcodes <- rownames(metadata_stored)[neighbours]

    # Compute proportion of neighbours in the barcode set
    sum(neighbour_barcodes %in% barcode_set) / kN
  }, FUN.VALUE = numeric(1))
}

.compile_pANN <- function(
    metadata_stored, neighbour_indices, ranges, kN, verbose
) {
  # Isolate columns & ensure cell names are present
  metadata_columns <- c("features", "counts", "mt.prop", "rb.prop",
                        "malat1", "Damage_level")

  metadata_plot <- metadata_stored[, metadata_columns]
  metadata_plot$Cells <- rownames(metadata_stored)
  damage_numbers <- as.numeric(unique(metadata_plot$Damage_level))
  damage_numbers <- damage_numbers[damage_numbers != 0]

  # Define sets of cells based on damage level
  barcodes <- list()
  for (number in damage_numbers) {
    barcode_name <- paste0("barcode_", number)
    barcodes[[barcode_name]] <- metadata_plot$Cells[
      metadata_plot$Damage_level == number
    ]
  }

  # Compute pANN for different damage levels
  if (verbose) {
    message("Computing pANN...")
  }

  # Compute proportion of artificial neighbours
  pANN_names <- c()
  for (number in damage_numbers) {
    barcode_name <- paste0("barcode_", number)
    pANN_name <- paste0("pANN_", number)
    pANN_names <- append(pANN_name, pANN_names)
    metadata_plot[[pANN_name]] <- .compute_pANN(
      barcodes[[barcode_name]], metadata_stored, neighbour_indices, kN
    )
  }

  # Find which damage population a cell is most frequently neighbouring
  # Isolate true cells & define relevant columns
  metadata_plot <- metadata_plot %>%
    dplyr::filter(.data$Damage_level == 0)

  # Find the column name with the maximum value for each true cell
  metadata_plot$max_pANN_col <- apply(
    metadata_plot[pANN_names], 1, function(row) pANN_names[which.max(row)]
  )

  return(metadata_plot)
}

.scale_damage <- function(metadata_plot, ranges, filter_threshold) {
  # Calculate min/max pANN value for each damage level
  pANN_summary_stats <- metadata_plot %>%
    dplyr::group_by(.data$max_pANN_col) %>%
    dplyr::summarise(
      min_value = min(get(first(.data$max_pANN_col)), na.rm = TRUE),
      max_value = max(get(first(.data$max_pANN_col)), na.rm = TRUE)
    )

  # Use the min/max values to scale pANN of each cell to lie between range
  metadata_plot <- metadata_plot %>%
    dplyr::mutate(
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
          expr(
            .data$max_pANN_col == !!name ~
              pANN_summary_stats$min_value[
                pANN_summary_stats$max_pANN_col == !!name
              ]
          )
        })
      ),
      max = case_when(
        !!!lapply(names(ranges), function(name) {
          expr(
            .data$max_pANN_col == !!name ~
              pANN_summary_stats$max_value[
                pANN_summary_stats$max_pANN_col == !!name
              ]
          )
        })
      ),
      value = case_when(
        !!!lapply(names(ranges), function(name) {
          expr(.data$max_pANN_col == !!name ~ get(!!name))
        })
      )
    )

  # Apply scaling
  metadata_plot$DamageDetective <- metadata_plot$lower_scale +
    ((metadata_plot$value - metadata_plot$min) /
       (metadata_plot$max - metadata_plot$min)) *
    (metadata_plot$upper_scale - metadata_plot$lower_scale)

  # Apply filter threshold to classify cells
  metadata_plot$DamageDetective_filter <- ifelse(
    metadata_plot$DamageDetective > filter_threshold, "damaged", "cell"
  )

  return(metadata_plot)
}

.finalise_output <- function(metadata_plot, count_matrix, filter_counts,
                             generate_plot, display_plot, target_damage,
                             filter_threshold, palette
) {
  # Apply filtering based on the DamageDetective score
  metadata_plot$DamageDetective_filter <- ifelse(
    metadata_plot$DamageDetective >= filter_threshold, "damaged", "cell"
  )
  output_columns <- c("Cells", "DamageDetective")
  metadata_output <- metadata_plot[, output_columns]

  # Filter the count matrix if requested
  if (filter_counts) {
    # Filter cells classified as "cell" (not damaged)
    metadata_filtered <- metadata_plot %>%
      dplyr::filter(.data$DamageDetective_filter == "cell")
    final_filtered_cells <- metadata_filtered$Cells
    final_filtered_matrix <- count_matrix[, final_filtered_cells]
    output <- final_filtered_matrix
  } else {
    # Only return damage labels
    output <- metadata_output
  }

  # Return both the output and the plot if requrested
  if (generate_plot) {
    final_plot <- plot_detection_outcome(
      qc_summary = metadata_plot,
      target_damage = target_damage,
      palette = palette
    )

    if (display_plot){
      print(final_plot)
    }

    output <- list(
      output = output,
      plot = final_plot
    )
  }
  return(output)
}
