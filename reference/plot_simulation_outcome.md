# plot_simulation_outcome

Function to generate quality control plots for altered and unaltered
counts.

## Usage

``` r
plot_simulation_outcome(
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

A `ggplot2` object representing the combined plot of altered and
unaltered counts.

## Details

This function generates a combined plot showing the distributions of
various quality control metrics (such as mitochondrial and ribosomal
proportions) before and after damage simulation. It compares unaltered
counts against altered counts for a more comprehensive assessment.
