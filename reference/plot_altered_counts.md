# plot_altered_counts

Function to generate a scatter plot of quality control metrics for
altered data.

## Usage

``` r
plot_altered_counts(
  qc_summary,
  palette = c("grey", "#7023FD", "#E60006"),
  target_damage = c(0.5, 1),
  base_size = 14
)
```

## Arguments

- qc_summary:

  A data frame containing the quality control summary for cells.

- palette:

  A character vector specifying the color gradient used for coloring the
  damage levels.

- target_damage:

  Numeric vector specifying the target damage levels for color scaling.

## Value

A list containing the plot object

## Details

This function visualizes the distribution of features and proportions of
mitochondrial and ribosomal genes for cells with altered counts. It also
adds a reference point to help assess the quality of altered data
against the expected distributions.
