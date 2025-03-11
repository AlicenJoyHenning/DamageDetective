#' detect_damage
#
#' Estimate the level of damage experienced by cells to inform the filtering
#' of cells that are most likely damaged where the amount of cytoplasmic RNA
#' lost is used as a proxy for damage. through comparison
#' to cells of the input data where RNA loss, ranging from 0 to 100 %, has
#' been simulated. The true and simulated cells are merged and processed
#' before the quality control metrics,
#'
#' * log(non-zero features)
#' * log(total counts)
#' * mitochondrial proportion
#' * ribosomal proportion
#' * MALAT1 expression
#'
#' are computed and compared through principal component analysis (PCA) to
#' generate a distance matrix. The top related cells, or nearest neighbours
#' (NNs), defined by the matrix are retrieved for each cell. For true cells,
#' the proportion of NNs that are artificial (pANNs), i.e. simulated cells,
#' are found for each level of loss simulated.The level of loss where the pANNs
#' is highest is used to assign a predicted level of loss to each true cell.
#' Damage labels are then assigned to true cells if the predicted level of
#' loss is greater than or equal to 40 % or a user specified threshold.
#'
#'
#' @param count_matrix Matrix or dgCMatrix containing the counts from
#'  single cell RNA sequencing data.
#' @param ribosome_penalty Numeric specifying the factor by which the
#'  probability of loosing a transcript from a ribosomal gene is multiplied by.
#'
#'  Values closer to 0 represent a greater penalty. We recommend using values
#'  around the suggested default as being more permissive of ribosomal loss
#'  leads to extensive reduction in ribosomal proportions and being overly
#'  restrictive of ribosomal loss leads to ribosomal proportions being
#'  relatively unchanged or even increased by the damage simulation, neither
#'  of which are observed in real single cell data.
#'
#'  The ideal penalty for your data can be found using the `select_penalty`
#'  function.
#'
#'  * Default is 0.01.
#' @param organism String specifying the organism of origin of the input
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
#' @param project_name String containing an identifier of the data being
#'  investigated to be used in the output plot titles.
#' @param filter_threshold Numeric between 0 and 1 specifying above what
#'   proportion of estimated cytoplasmic RNA loss a cell should be regarded
#'   as damaged. It is recommended not going lower than 0.4 as there is often no
#'   considerable difference to true cells when less than 0.4 of the total
#'   cytoplasmic RNA is lost.
#'
#'  * Default 0.75.
#' @param damage_levels Specification of how many sets of artificial damaged
#'  cells are created. The true cells will be compared to cells from each
#'  set and the set with which it has the greatest proportion of nearest
#'  neighbours will become the estimated level of damage for the true cell.
#'
#'  There are 3 default options provided,
#'  - 3
#'  - 5
#'  - 7
#'
#' Which include the following sets of damaged cells,
#' - 3: c(0.00001, 0.08), c(0.1, 0.4), c(0.5, 0.9)
#' - 5: c(0.00001, 0.08), c(0.1, 0.3), c(0.3, 0.5), c(0.5, 0.7), c(0.7, 0.9)
#' - 7: c(0.00001, 0.08), c(0.1, 0.3), c(0.3, 0.4), c(0.4, 0.5), c(0.5, 0.7),
#' c(0.7, 0.9), c(0.9, 0.99999)
#'
#' Else, a user is also able to specify their own sets as follows, though
#' we highly recommend the use of the defaults.
#'     damage_levels = list(
#'     pANN_50 = c(0.01, 0.5),
#'     pANN_90 = c(0.5, 0.9)
#' )
#'
#' * Default is 5.
#' @param mito_quantile Numeric specifying below what proportion of
#'  mitochondrial content cells are used for sampling for simulation.
#'
#' * Default is 0.75, meaning only cells with less than 0.75 proportion of
#' mitochondrial counts are sampled for simulated.
#' @param kN Numeric specifying the number of neighbours counted when
#'  calculating the proportion of artificial nearest neighbours of each
#'  true cell. It is recommended for this value to correspond to
#'  the number of cells present in your data.
#'
#' * Default is a third of the number of cells present in the input matrix.
#' @param include_pANN Boolean specifying whether the function output
#'    should include the proportion of nearest neighbours for every cell
#'    to each set of artificially damaged cells. This can be used to
#'    override the selection of the most suitable set of artificial cells
#'    to describe each true cell.
#'
#'  * Default is FALSE.
#' @param filter_counts Boolean specifying whether to output only the true
#'  cells whose proportion of artificial nearest neighbours showed them
#'  to be most similar to artificially damaged cells that had lost less RNA
#'  than that specified in the filter_threshold. This automatically filters
#'  the data for a user.
#'
#'  If FALSE, the output provided is a data frame containing the barcodes
#'  of the input data alongside their estimated levels of damage and a
#'  classification as either 'damaged' or 'cell' based on the filter_threshold.
#'  This allows a user control over filtering decisions.
#'
#'  * Default is FALSE.
#' @param verbose Boolean specifying whether messages and function progress
#'  should be displayed in the console.
#'
#'  * Default is TRUE.
#' @return Either a filtered count matrix or a data frame with annotations
#'  of estimated damage level for each cell and a classification as either
#'  'damaged' or 'cell'.
#' @export
#' @examples
#' data("test_counts", package = "DamageDetective")
#'
#' test <- detect_damage(
#'   count_matrix = test_counts,
#'   ribosome_penalty = 0.05,
#'   filter_threshold = 0.75,
#'   damage_levels = 3,
#'   kN = 100,
#'    project_name = "Test data",
#'    save_plot = temp_dir*()
#'    )
#' @importFrom dplyr %>% group_by summarise mutate case_when
#' @importFrom ggplot2 ggsave theme element_rect margin
#' @importFrom cowplot ggdraw draw_label plot_grid
#' @importFrom ggpubr get_legend
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom Seurat ScaleData PercentageFeatureSet FetchData
#' @importFrom rlang expr !!!
#' @importFrom stats prcomp
detect_damage <- function(
    count_matrix,
    ribosome_penalty,
    organism = "Hsap",
    project_name = "Project_name",
    filter_threshold = 0.5,
    damage_levels = 5,
    mito_quantile = 0.75,
    kN = NULL,
    include_pANN = FALSE,
    filter_counts = FALSE,
    verbose = TRUE,
    ...
){
  # Phase One: data preparations for simulation ----

  # Retrieve genes corresponding to the organism of interest
  # Human
  if (organism == "Hsap") {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^(RPS|RPL)"
    nuclear <- c("FIRRE", "NEAT1","XIST", "MALAT1", "MIAT", "MEG3", "KCNQ1OT1", "HOXA11-AS", "FTX")
    MALAT1 <- "MALAT1"
  }

  # Mouse
  if (organism == "Mmus") {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^(rps|rpl)"
    nuclear <- c("Firre", "Neat1","Xist", "Malat1", "Miat", "Meg3", "Kcnq1ot1", "Hoxa11-as", "Ftx")
    MALAT1 <- "Malat1"
  }

  # Allow for user specification for non-standard organism
  if (!organism %in% c("Hsap", "Mmus")) {
    mito_pattern <- organism$mito_pattern
    ribo_pattern <- organism$ribo_pattern
    nuclear <- organism$nuclear
  }

  # Collect mito & ribo information of all true cells
  mito_idx <- grep(mito_pattern, rownames(count_matrix), ignore.case = FALSE)
  ribo_idx <- grep(ribo_pattern, rownames(count_matrix), ignore.case = FALSE)
  mito <- colSums(count_matrix[mito_idx, , drop = FALSE]) / colSums(count_matrix)
  ribo <- colSums(count_matrix[ribo_idx, , drop = FALSE]) / colSums(count_matrix)

  df <- data.frame(mito = mito,
                   ribo = ribo,
                   row.names = colnames(count_matrix))

  # Identify high mito proportion cells (likely damaged)
  mito_top <- quantile(df$mito, probs = mito_quantile, na.rm = TRUE)

  # Subset true cells below the 75th percentile
  filtered_cells <- rownames(df)[df$mito <= mito_top]
  matrix_filtered <- count_matrix[, filtered_cells]

  # Phase Two: Generate simulations for defined damage levels ----

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
    level <- names(ranges)[i]  # e.g., "pANN_10"
    damage_range <- ranges[[i]]
    damage_ranges <- append(ranges[[i]], damage_ranges)

    if (verbose){
      message("Simulating data cells between ", damage_range[[1]], " and ", damage_range[[2]], " level of damage...")
    }

    # Simulate counts
    damaged_cells <- simulate_counts(
      count_matrix = matrix_filtered,
      damage_proportion = 0.15,
      target_damage = damage_range,
      ribosome_penalty = ribosome_penalty,
      damage_distribution = "symmetric",
      damage_steepness = "steep",
      save_plot = NULL
    )

    # Isolate counts of artificial cells (newly damaged)
    barcodes <- subset(damaged_cells$qc_summary, Damaged_Level != 0)$Cell
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
  combined$mt.prop <- PercentageFeatureSet(combined, pattern = mito_pattern) / 100
  combined$rb.prop <- PercentageFeatureSet(combined, pattern = ribo_pattern) / 100
  combined$malat1 <- FetchData(combined, vars = MALAT1)

  # Phase Three: Perform PCA on the meta data values ----

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
  metadata <- subset(metadata, Damage_level == 0)

  # Find the column name with the maximum value for each true cell
  metadata$max_pANN_col <- apply(metadata[pANN_names], 1, function(row) {
    pANN_names[which.max(row)]
  })

  # Calculate min/max pANN value for each damage level
  pANN_summary_stats <- metadata %>%
    dplyr::group_by(max_pANN_col) %>%
    dplyr::summarise(
      min_value = min(get(first(max_pANN_col)), na.rm = TRUE),
      max_value = max(get(first(max_pANN_col)), na.rm = TRUE)
    )

  # Use the min/max values to scale to pANN of each cell to lie between range
  metadata <- metadata %>%
    mutate(
      lower_scale = case_when(
        !!!lapply(names(ranges), function(name) {
          rlang::expr(max_pANN_col == !!name ~ ranges[[!!name]][[1]])
        })
      ),
      upper_scale = case_when(
        !!!lapply(names(ranges), function(name) {
          rlang::expr(max_pANN_col == !!name ~ ranges[[!!name]][[2]])
        })
      ),
      min = case_when(
        !!!lapply(names(ranges), function(name) {
          rlang::expr(max_pANN_col == !!name ~ pANN_summary_stats$min_value[pANN_summary_stats$max_pANN_col == !!name])
        })
      ),
      max = case_when(
        !!!lapply(names(ranges), function(name) {
          rlang::expr(max_pANN_col == !!name ~ pANN_summary_stats$max_value[pANN_summary_stats$max_pANN_col == !!name])
        })
      ),
      value = case_when(
        !!!lapply(names(ranges), function(name) {
          rlang::expr(max_pANN_col == !!name ~ get(!!name))
        })
      )
    )

  # Apply scaling
  metadata$scaled_pANN <- metadata$lower_scale + ((metadata$value  - metadata$min) / (metadata$max - metadata$min)) * (metadata$upper_scale - metadata$lower_scale)


  # Phase Four: Filter, mark and visualise damaged cells  ----
  metadata$DamageDetective <- ifelse(metadata$scaled_pANN >= filter_threshold, "damaged", "cell")

  if (include_pANN){
    # Include the proportions for each damage level simulated
    columns <- c("Cells", "scaled_pANN", "DamageDetective")
    columns <- append(pANN_names, columns)
    metadata_output <- metadata[, columns]
  } else {
    # Return only the scaled proportions
    metadata_output <- metadata[, c("Cells", "scaled_pANN", "DamageDetective")]
  }

  # Visualise cells according to damage level
  if (!is.null(save_plot)){

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

    # Save the final plot to a specified directory
    output_file <- file.path(save_plot, paste0("DamageDetective.", plot_type))

    # Use ggsave to save the plot
    if (plot_type == "png") {
      ggplot2::ggsave(
        output_file,
        plot = final_plot,
        width = 8,
        height = 5,
        dpi = 300,
        units = "in"
      )
    } else if (plot_type == "svg") {
      ggplot2::ggsave(
        output_file,
        plot = final_plot,
        width = 8,
        height = 5,
        units = "in"
      )
    }

    output <- list(
      metadata = metadata_output,
      plot = final_plot
    )

  }

  # If specified, filter and return the count matrix only
  if (filter_counts){
    metadata_filtered <- subset(metadata, DamageDetective == "cell")
    final_filtered_cells <- metadata_filtered$Cells
    final_filtered_matrix <- count_matrix[, final_filtered_cells]
    output <- final_filtered_matrix
  }

  return(output)

}
