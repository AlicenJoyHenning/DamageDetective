
#' detect_damage
#
#' Estimate the level of cytoplasmic RNA loss in a cell, a fundamental
#' proxy for damage in single cell RNA sequencing, through comparison
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
#' @param filter_threshold Numeric between 0 and 1 specifying above what
#'   proportion of estimated cytoplasmic RNA loss a cell should be regarded
#'   as damaged. We recommend not going lower than 0.4.
#'
#'  * Default 0.4
#' @param verbose Boolean specifying whether messages and function progress
#'  should be displayed in the console.
#'
#'  * Default is TRUE.
#' @return
#' @export
#'
#' @examples

# 1. Create df with mito, ribo information of all true cells
mito_idx <- grep("^MT-", rownames(matrix), ignore.case = FALSE)
ribo_idx <- grep("^(RPS|RPL)", rownames(matrix), ignore.case = FALSE)
mito <- colSums(matrix[mito_idx, , drop = FALSE]) / colSums(matrix)
ribo <- colSums(matrix[ribo_idx, , drop = FALSE]) / colSums(matrix)

df <- data.frame(mito = mito,
                 ribo = ribo,
                 row.names = colnames(matrix))

# Identify high mito proportion cells (likely damaged)
mito_quantile = 0.75
mito_top <- quantile(df$mito, probs = mito_quantile, na.rm = TRUE)

# 2. Subset to true cells below the 75th percentile
filtered_cells <- rownames(df)[df$mito <= mito_top]
matrix <- matrix[, filtered_cells]

# Define column-specific ranges
ranges <- list(
  pANN_10 = c(0.00001, 0.08),
  pANN_30 = c(0.1, 0.3),
  pANN_50 = c(0.3, 0.5),
  pANN_70 = c(0.5, 0.7),
  pANN_90 = c(0.7, 0.9)
)

# Define corresponding damage distributions
damage_distributions <- c("symmetric", "left_skewed", "symmetric", "symmetric", "left_skewed")

# Initialize matrix with original counts
matrix_updated <- matrix

# Iterate through ranges and simulate damage
for (i in seq_along(ranges)) {
  level <- names(ranges)[i]  # e.g., "pANN_10"
  damage_range <- ranges[[i]]
  damage_dist <- damage_distributions[i]

  # Simulate counts
  damaged_cells <- simulate_counts(
    count_matrix = matrix,
    damage_proportion = 0.15,
    target_damage = damage_range,
    ribosome_penalty = penalty,
    damage_distribution = damage_dist,
    damage_steepness = "steep"
  )

  # Isolate counts of damaged cells
  barcodes <- subset(damaged_cells$qc_summary, Damaged_Level != 0)$Cell
  damaged_matrix <- damaged_cells$matrix[, barcodes]

  # Append the level of damage to the barcodes
  colnames(damaged_matrix) <- paste0(colnames(damaged_matrix), "_", gsub("pANN_", "", level))

  # Merge with updated matrix
  matrix_updated <- cbind(matrix_updated, damaged_matrix)

}


# Combine with true cells and remove duplicated, unaltered true cells
matrix_combined <- cbind(matrix_updated, matrix) # must be original matrix (without empty droplets)
matrix_combined <- matrix_combined[, unique(colnames(matrix_combined))]

# Combine artificial cells with true cells ------

# Create Seurat object with all cells
combined <- CreateSeuratObject(counts = matrix_combined)

# Pre-processing the merged object
combined <- NormalizeData(combined) %>%
  FindVariableFeatures() %>%
  ScaleData()

# Add meta_data label with level of damage of the cell
combined$Damage_level <- sub(".*_", "", rownames(combined@meta.data))
combined$Damage_level <- ifelse(!combined$Damage_level %in% c("10", "30", "50", "70", "90"), 0, combined$Damage_level)
# table(combined$Damage_level) # Verify different levels

# Find mitochondrial and ribosomal percentage for each cell
combined$log.features <- log10(combined$nFeature_RNA)
combined$log.counts <- log10(combined$nCount_RNA)
combined$mt.prop <- PercentageFeatureSet(combined, pattern = "^MT-") / 100
combined$rb.prop <- PercentageFeatureSet(combined, pattern = "^(RPS|RPL)") / 100
combined$malat1 <- FetchData(combined, vars = "MALAT1")


# PCA on the meta data values ----

# Isolate variables of interest
metadata <- combined@meta.data
metadata <- metadata[, c("log.features", "log.counts", "mt.prop", "rb.prop", "malat1")]
# metadata <- metadata[, c("nFeature_RNA", "nCount_RNA", "mt.prop", "rb.prop", "malat1")]

# Perform PCA
pca_result <- prcomp(metadata, center = TRUE, scale. = TRUE)

# Extract PCA embeddings of top principal components
pca_coord <- pca_result$x[, 1:5]

# Calculate the euclidean distance between PC embeddings of cells
dist_mat <- rdist(pca_coord)  # Compute Euclidean distances
rownames(dist_mat) <- rownames(metadata)  # Assign row names
colnames(dist_mat) <- rownames(metadata)  # Assign column names

# Isolate columns of for PCA & ensure cell names present
#metadata <- combined@meta.data[, c("log.features", "log.counts", "mt.prop", "rb.prop", "malat1", "Damage_level")]
metadata <- combined@meta.data[, c("nFeature_RNA", "nCount_RNA", "mt.prop", "rb.prop", "malat1", "Damage_level")]

metadata$Cells <- rownames(metadata)

# Define sets of cells based on damage level
barcodes_10 <- metadata$Cells[metadata$Damage_level == "10"]
barcodes_30 <- metadata$Cells[metadata$Damage_level == "30"]
barcodes_50 <- metadata$Cells[metadata$Damage_level == "50"]
barcodes_70 <- metadata$Cells[metadata$Damage_level == "70"]
barcodes_90 <- metadata$Cells[metadata$Damage_level == "90"]

# Find the proportion of nearest neighbours belonging to a damage level
# Function to compute pANN for a given barcode set
compute_pANN <- function(barcode_set) {
  sapply(rownames(dist_mat), function(cell) {
    # Find indices of the 1000 nearest neighbors (excluding itself)
    neighbors <- order(dist_mat[cell, ])[2:1001]
    neighbor_barcodes <- rownames(dist_mat)[neighbors]

    # Compute proportion of neighbors in the barcode set
    sum(neighbor_barcodes %in% barcode_set) / 1000
  })
}

# Compute pANN for different damage levels
metadata$pANN_10 <-compute_pANN(barcodes_10)
metadata$pANN_30 <- compute_pANN(barcodes_30)
metadata$pANN_50 <- compute_pANN(barcodes_50)
metadata$pANN_70 <- compute_pANN(barcodes_70)
metadata$pANN_90 <- compute_pANN(barcodes_90)

# Find which damage population a cell is most frequently neighboring
# Isolate  true cells & define relevant columns
metadata <- subset(metadata, Damage_level == 0)
pANN_cols <- c("pANN_10", "pANN_30", "pANN_50", "pANN_70", "pANN_90")

# Find the column name with the maximum value for each true cell
metadata$max_pANN_col <- apply(metadata[pANN_cols], 1, function(row) {
  pANN_cols[which.max(row)]
})

# Calculate min/max for the cells belonging to each pANN group
pANN_summary_stats <- metadata %>%
  dplyr::group_by(max_pANN_col) %>%
  dplyr::summarise(
    min_value = min(get(first(max_pANN_col)), na.rm = TRUE),
    max_value = max(get(first(max_pANN_col)), na.rm = TRUE)
  )


metadata <- metadata %>%
  mutate(
    lower_scale = case_when(
      max_pANN_col == "pANN_10" ~ ranges$pANN_10[[1]],
      max_pANN_col == "pANN_30" ~ ranges$pANN_30[[1]],
      max_pANN_col == "pANN_50" ~ ranges$pANN_50[[1]],
      max_pANN_col == "pANN_70" ~ ranges$pANN_70[[1]],
      max_pANN_col == "pANN_90" ~ ranges$pANN_90[[1]]
    ),
    upper_scale = case_when(
      max_pANN_col == "pANN_10" ~ ranges$pANN_10[[2]],
      max_pANN_col == "pANN_30" ~ ranges$pANN_30[[2]],
      max_pANN_col == "pANN_50" ~ ranges$pANN_50[[2]],
      max_pANN_col == "pANN_70" ~ ranges$pANN_70[[2]],
      max_pANN_col == "pANN_90" ~ ranges$pANN_90[[2]]
    ),
    min = case_when(
      max_pANN_col == "pANN_10" ~ pANN_summary_stats$min_value[pANN_summary_stats$max_pANN_col == "pANN_10"],
      max_pANN_col == "pANN_30" ~ pANN_summary_stats$min_value[pANN_summary_stats$max_pANN_col == "pANN_30"],
      max_pANN_col == "pANN_50" ~ pANN_summary_stats$min_value[pANN_summary_stats$max_pANN_col == "pANN_50"],
      max_pANN_col == "pANN_70" ~ pANN_summary_stats$min_value[pANN_summary_stats$max_pANN_col == "pANN_70"],
      max_pANN_col == "pANN_90" ~ pANN_summary_stats$min_value[pANN_summary_stats$max_pANN_col == "pANN_90"]
    ),
    max = case_when(
      max_pANN_col == "pANN_10" ~ pANN_summary_stats$max_value[pANN_summary_stats$max_pANN_col == "pANN_10"],
      max_pANN_col == "pANN_30" ~ pANN_summary_stats$max_value[pANN_summary_stats$max_pANN_col == "pANN_30"],
      max_pANN_col == "pANN_50" ~ pANN_summary_stats$max_value[pANN_summary_stats$max_pANN_col == "pANN_50"],
      max_pANN_col == "pANN_70" ~ pANN_summary_stats$max_value[pANN_summary_stats$max_pANN_col == "pANN_70"],
      max_pANN_col == "pANN_90" ~ pANN_summary_stats$max_value[pANN_summary_stats$max_pANN_col == "pANN_90"]
    ),
    value = case_when(
      max_pANN_col == "pANN_10" ~ pANN_10,
      max_pANN_col == "pANN_30" ~ pANN_30,
      max_pANN_col == "pANN_50" ~ pANN_50,
      max_pANN_col == "pANN_70" ~ pANN_70,
      max_pANN_col == "pANN_90" ~ pANN_90
    )
  )

# Apply scaling
metadata$scaled_pANN <- metadata$lower_scale + ((metadata$value  - metadata$min) / (metadata$max - metadata$min)) * (metadata$upper_scale - metadata$lower_scale)

# Visualise cells according to damage level ----

metadata$scaled_pANN <- metadata$scaled_pANN # / 100 # Scale for colour gradient

# Looking at the level of isolated subsets to see if things have gone well
high_subset <- subset(metadata, scaled_pANN >= 0.5)
low_subset <- subset(metadata, scaled_pANN < 0.5)

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

plot_mito_ribo | plot_mito_features
