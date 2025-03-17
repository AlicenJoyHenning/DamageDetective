# Function to randomly remove transcripts
  apply_transcript_loss <- function(
    matrix, 
    cell_indices = damaged_cell_selections, 
    gene_indices = non_mito_idx, 
    total_loss = damage_label$damage_level) {
    
    # Loop over each cell
    for (cell in cell_indices) {
      
      # Skip cells where no loss is required
      loss_proportion <- total_loss[cell]
      loss_target <- round(loss_proportion * sum(matrix[, cell]))

      if (loss_target == 0) next
      
      # Extract counts for target genes in this cell
      gene_counts <- matrix[gene_indices, cell]

      # Compute probability of loss for each gene (proportional to count)
      probabilities <- gene_counts / sum(gene_counts)
      
      # Sample transcript losses based on probabilities
      loss_sample <- sample(gene_indices, size = loss_target, replace = TRUE, prob = probabilities)
      
      # Collect the total transcripts to lose from each gene
      loss_counts <- table(loss_sample)
      
      # Update the matrix by subtracting loss counts from the corresponding genes
      for (gene in names(loss_counts)) {
        gene_index <- as.integer(gene)  
        matrix[gene_index, cell] <- pmax(0, matrix[gene_index, cell] - loss_counts[gene])
      }
    }
    
    return(matrix)
  }
  
  
  # Apply function
  damaged_count_matrix <- apply_transcript_loss(damaged_count_matrix)
  