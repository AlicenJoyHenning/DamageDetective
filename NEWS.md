# DamageDetective (development version)

# DamageDetective 0.99.0

- The `DamageDetective` package has been released and is available for  
  installation from GitHub. CRAN submission is in progress.

# DamageDetective 1.0.0  

- CRAN submission comments have been received and addressed.  
- Fixed DOI formatting.  
- Fixed spacing in the vignette.  
- Ensured code lines are shorter than 80 characters for good practice.  

# DamageDetective 2.0.0  

- Replaced the `.perturb_cells()` function with a C++ implementation to preserve transcript-level sampling while significantly improving runtime. Core functionality remains the same, with the exception of C++-based indexing (0-based vs. 1-based in R).
- Introduced a clustering-based damage detection approach: cells are grouped into populations using low-resolution Seurat clustering, and extreme cytoplasmic RNA loss is simulated within each group to better model realistic damaged profiles.
- Updated documentation throughout the package to reflect changes in methodology, implementation, and usage.

# DamageDetective 2.0.1
- New feature: DamageDetective now checks entire clusters for evidence of damage before simulating artificial damage. Clusters are flagged as damaged if their top differentially expressed genes are predominantly mitochondrial. This suggests a lack of true biological signal, a hallmark of heavily damaged cell populations.

- Damage score assignment: As a result, cells are now assigned a damage score of 1 if they belong to a cluster identified as likely damaged. Cells not flagged in this way continue to receive a score based on proximity to simulated damaged cells in PC space.
