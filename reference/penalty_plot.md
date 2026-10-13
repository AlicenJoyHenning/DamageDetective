# Penalty Plot

A precomputed `ggplot` object visualizing the effects of different
ribosomal penalties on simulated damaged cell populations.

## Usage

``` r
data(penalty_plot)
```

## Format

A `gg` object created with the `ggplot2` and `patchwork` packages.

## Details

This object consists of a combined patchwork plot with three panels: one
for each ribosomal penalty level (no penalty, medium penalty, heavy
penalty), generated using
[`DamageDetective::simulate_counts()`](https://alicenjoyhenning.github.io/DamageDetective/reference/simulate_counts.md).
These plots were created to illustrate how ribosomal content
penalization influences the separation between healthy and damaged
cells.

## See also

[`simulate_counts`](https://alicenjoyhenning.github.io/DamageDetective/reference/simulate_counts.md)
