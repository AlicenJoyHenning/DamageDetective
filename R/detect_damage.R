# Function to perf# Data preparation ----
# Load the matrix without cell types labelled 
PH_CT6_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_high_CT6_reference_matrix.rds")
colnames(PH_CT6_matrix) <- names(colnames(PH_CT6_matrix))

# Running simulation ----
# Generate artificial datasets covering range of damage levels 

# Filter the input data to contain cells where mito.prop is not far from mean 
mito_idx <- grep("^MT-", rownames(PH_CT6_matrix), ignore.case = FALSE)
mito <- colSums(PH_CT6_matrix[mito_idx, ]) / colSums(PH_CT6_matrix)
IQR <- IQR(mito)
median <- median(mito)
filtered_cells <- names(mito[(mito < (median + IQR)) & (mito > (median - IQR))])
PH_CT6_filtered_matrix <- PH_CT6_matrix[, filtered_cells]

# Define column-specific ranges
ranges <- list(
  pANN_10 = c(0.00001, 0.08),
  pANN_30 = c(0.1, 0.3),
  pANN_50 = c(0.3, 0.5),
  pANN_70 = c(0.5, 0.7),
  pANN_90 = c(0.7, 0.9)
)


# < 10 %
damaged_cells_10 <- simulate_counts(
  count_matrix = PH_CT6_filtered_matrix,
  damage_proportion = 0.15,
  target_damage = ranges$pANN_10,
  damage_distribution = "symmetric",
  damage_steepness = "steep"
)

# Isolate counts of damaged cells
barcodes_10 <- subset(damaged_cells_10$qc_summary, Damaged_Level != 0)
barcodes_10 <- barcodes_10[, c("Cell")]
damaged_10_matrix <- damaged_cells_10$matrix[, barcodes_10]
# Append the level of damage to the barcodes of the damaged cells
colnames(damaged_10_matrix) <- paste0(colnames(damaged_10_matrix), "_10")

# Merge the counts with the original
PH_CT6_matrix_updated <- cbind(PH_CT6_matrix, damaged_10_matrix)


# 10 < 30 %
damaged_cells_30 <- simulate_counts(
  count_matrix = PH_CT6_filtered_matrix, 
  damage_proportion = 0.15,
  target_damage = ranges$pANN_30,
  damage_distribution = "left_skewed",
  damage_steepness = "steep"
)

# Isolate counts of damaged cells 
barcodes_30 <- subset(damaged_cells_30$qc_summary, Damaged_Level != 0)
barcodes_30 <- barcodes_30[, c("Cell")]
damaged_30_matrix <- damaged_cells_30$matrix[, barcodes_30]
# Append the level of damage to the barcodes of the damaged cells
colnames(damaged_30_matrix) <- paste0(colnames(damaged_30_matrix), "_30")

# Merge the counts with the original 
PH_CT6_matrix_updated <- cbind(PH_CT6_matrix_updated , damaged_30_matrix)

# 30 < 50 %
damaged_cells_50 <- simulate_counts(
  count_matrix = PH_CT6_filtered_matrix, 
  damage_proportion = 0.15,
  target_damage = ranges$pANN_50,
  damage_distribution = "symmetric",
  damage_steepness = "steep"
)

# Isolate counts of damaged cells 
barcodes_50 <- subset(damaged_cells_50$qc_summary, Damaged_Level != 0)
barcodes_50 <- barcodes_50[, c("Cell")]
damaged_50_matrix <- damaged_cells_50$matrix[, barcodes_50]
# Append the level of damage to the barcodes of the damaged cells
colnames(damaged_50_matrix) <- paste0(colnames(damaged_50_matrix), "_50")

# Merge the counts with the original 
PH_CT6_matrix_updated <- cbind(PH_CT6_matrix_updated, damaged_50_matrix)

# 50 < 70 %
damaged_cells_70 <- simulate_counts(
  count_matrix = PH_CT6_filtered_matrix, 
  damage_proportion = 0.15,
  target_damage = ranges$pANN_70,
  damage_distribution = "symmetric",
  damage_steepness = "steep"
)

# Isolate counts of damaged cells 
barcodes_70 <- subset(damaged_cells_70$qc_summary, Damaged_Level != 0)
barcodes_70 <- barcodes_70[, c("Cell")]
damaged_70_matrix <- damaged_cells_70$matrix[, barcodes_70]
# Append the level of damage to the barcodes of the damaged cells
colnames(damaged_70_matrix) <- paste0(colnames(damaged_70_matrix), "_70")

# Merge the counts with the original 
PH_CT6_matrix_updated <- cbind(PH_CT6_matrix_updated, damaged_70_matrix)

# 70 < 90 %
damaged_cells_90 <- simulate_counts(
  count_matrix = PH_CT6_filtered_matrix, 
  damage_proportion = 0.15,
  target_damage = ranges$pANN_90,
  damage_distribution = "symmetric", # Too high up to 100 % not as likely
  damage_steepness = "steep"
)

# Isolate counts of damaged cells 
barcodes_90 <- subset(damaged_cells_90$qc_summary, Damaged_Level != 0)
barcodes_90 <- barcodes_90[, c("Cell")]
damaged_90_matrix <- damaged_cells_90$matrix[, barcodes_90]
# Append the level of damage to the barcodes of the damaged cells
colnames(damaged_90_matrix) <- paste0(colnames(damaged_90_matrix), "_90")

# Merge the counts with the original 
PH_CT6_matrix_updated <- cbind(PH_CT6_matrix_updated, damaged_90_matrix)


# Combining artificial cells with true cells ------

# Create Seurat object with all cells
combined <- CreateSeuratObject(counts = PH_CT6_matrix_updated )

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

# Perform PCA 
pca_result <- prcomp(metadata, center = TRUE, scale. = TRUE)

# Extract PCA embeddings of top principal components 
pca_coord <- pca_result$x[, 1:5]

# Calculate the euclidean distance between PC embeddings of cells 
dist_mat <- rdist(pca_coord)  # Compute Euclidean distances
rownames(dist_mat) <- rownames(metadata)  # Assign row names
colnames(dist_mat) <- rownames(metadata)  # Assign column names

# Isolate columns of for PCA & ensure cell names present
metadata <- combined@meta.data[, c("log.features", "log.counts", "mt.prop", "rb.prop", "malat1", "Damage_level")]
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
  x = "log.features", 
  y = "mt.prop",
  damage_column = "scaled_pANN",
  altered = TRUE,
  target_damage = c(0.5, 0.9)
)


plot_mito_ribo | plot_mito_features


# Set threshold for filtering 
threshold <- 0.4
metadata$DamageDetective <- ifelse(metadata$scaled_pANN >= threshold, "damaged", "cell")
metadata_filtered <- subset(metadata, DamageDetective == "cell")


plot_mito_ribo <- plot_outcome(
  data = metadata_filtered, 
  x = "rb.prop", 
  y = "mt.prop",
  damage_column = "scaled_pANN",
  altered = TRUE,
  mito_ribo = TRUE,
  target_damage = c(0.5, 0.9)
)

plot_mito_features <- plot_outcome(
  data = metadata_filtered, 
  x = "log.features", 
  y = "mt.prop",
  damage_column = "scaled_pANN",
  altered = TRUE,
  target_damage = c(0.5, 0.9)
)


plot_mito_ribo | plot_mito_features
orm the damaged cell predictions
