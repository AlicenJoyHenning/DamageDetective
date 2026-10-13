# plot_unaltered_counts

Function to generate a scatter plot of quality control metrics for
unaltered data.

## Usage

``` r
plot_unaltered_counts(qc_summary, base_size = 14)
```

## Arguments

- qc_summary:

  A data frame containing the quality control summary for cells.

## Value

A `ggplot2` object representing the scatter plot of quality control
metrics for unaltered cells.

## Details

This function visualizes the distribution of features and proportions of
mitochondrial and ribosomal genes for unaltered cells. It provides a
reference plot to assess the original data's quality before any
alterations.
