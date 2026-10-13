# plot_ribosomal_penalty

Function to generate a scatter plot of the simulated data focusing on
ribosomal proportion.

## Usage

``` r
plot_ribosomal_penalty(
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

A ggplot2 plot

## Details

This function visualizes the distribution of mitochondrial and ribosomal
proportion for cells with altered counts.
