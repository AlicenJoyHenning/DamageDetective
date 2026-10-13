# plot_detection_outcome

Function to generate a plot showing the distribution of quality control
metrics across altered data.

## Usage

``` r
plot_detection_outcome(
  qc_summary,
  target_damage = c(0.1, 0.8),
  palette = c("grey", "#7023FD", "#E60006"),
  base_size = 14
)
```

## Arguments

- qc_summary:

  A data frame containing the quality control summary for cells.

- target_damage:

  Numeric vector specifying the target damage levels for color scaling.

- palette:

  A character vector specifying the color gradient used for coloring the
  damage levels.

## Value

A `ggplot2` object representing the scatter plot of quality control
metrics.

## Details

This function visualizes the distribution of features and proportions of
mitochondrial and ribosomal genes for altered cells, coloring the points
by their damage levels. It helps in assessing how well the damage
detection process has classified cells based on their quality control
metrics.
