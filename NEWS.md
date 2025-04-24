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
