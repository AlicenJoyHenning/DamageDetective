# DamageDetective <img src="man/figures/logo.svg" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cosimameyer/overviewR/actions) ![Build Status](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/build.yml/badge.svg)

<!-- badges: end -->

## Content

[Description](#description) \| [Installation](#installation) \| [Quick start](#quick-start) \| [Authors](#authors) \| [License](#license) \| [References](#references)

## Description

Damaged cells are a class of low-quality artifacts in single-cell RNA sequencing (scRNA-seq) describing cells that have succumbed to stress before being sequenced. As this results in altered gene expression profiles, retaining damaged cells in downstream analyses compromises the reliability of the analysis. The detection and removal of damaged cells is therefore a standard practice in the pre-processing of scRNA-seq data.

Current approaches detect damage according to deviations in cell-level quality control metrics influenced by the loss of plasma membrane integrity, a well-accepted consequence of damage. However, this approach assumes all viable cells follow similar distributions, an assumption that does not always hold in heterogeneous samples. More recent methods improve upon this by analysing cells at a population level, isolating cells with similar distributions before assessing finer deviations. However, this assumes that all distinct populations represent true cells meaning if abundant, damaged cells risk misclassification as true cells. Ultimately, the filtering decisions of current approaches are motivated more by statistical definitions of deviation than biological definitions of damage, a non-linear, stochastic process that varies across cell types. 

`DamageDetective` takes a different approach inspired by `DoubletFinder`$^*$, a high performing community accepted tool for doublet removal in scRNA-seq. Here, rather than detecting damage by measuring the extent to which cells deviate from each other, it is detected by measuring the extent to which cells align with artificially damaged versions of themselves. Grounded, too, by loss of plasma membrane integrity, this approach simulates damage through the probabilistic escape of cytoplasmic RNA where the proportion of escape is assumed to directly reflect damage severity. By comparing the expression profiles of true cells to artificially damaged cell profiles in reduced-dimensional space, `DamageDetective` estimates the damage severity of true cells in a score from 0 to 1. This provides an intuitive scale for filtering, based directly on biological definitions of damage, that is standardised across cell types, samples, and experiments.

## Installation

`DamageDetective` can be installed from CRAN using,

``` r
install.packages('DamageDetective')
```

Alternatively, the latest development version can be installed directly from GitHub using,

``` r
library(devtools)
devtools::install_github("AlicenJoyHenning/DamageDetective", build_vignettes = TRUE)
```

To verify installation, run the following to see if you can view the package vignette and the function help page,

``` r
library(DamageDetective)
help(package = "DamageDetective")
?DamageDetective()
```

------------------------------------------------------------------------

## Quick start

The demonstrations below can be followed immediately after loading the package and serve as a test to ensure all is running smoothly. For function descriptions and usage examples please refer to the [package vignette](link).

Begin by loading the dummy data provided by the package, `test_counts`, an artificial PBMC dataset containing 500 cells and 10009 genes.

``` r
data("test_counts", package = "DamageDetective")
dim(test_counts)
# [1] 10000   500
```

<ul>

<li>

<h3>Return default output</h3>

</li>

</ul>

The primary goal of `DamageDetective` is to inform the filtering of damaged cells from single cell data. This can be achieved using the `detect_damage` function that requires a count matrix as input and returns a data frame containing the barcodes of the count matrix with the estimated levels of damage.

The results can be accessed from the `output` slot and used to inform filtering of cells from the count matrix. Below we are filtering cells with a damage level above 70 %,

``` r
# Perform damage detection
default_test <- detect_damage(
  count_matrix = test_counts
)

# View output
head(default_test$qc_summary) 
# 

# Filter cells with an estimated damage level above 70 % 
undamaged_cells <- subset(default_test$qc_summary, DamageDetective > 0.7)
filtered_counts <- test_counts[, undamaged_cells]
```

> By default, `detect_damage` will provide plots of the data where each cell is a point coloured according to the damage level estimated by `DamageDetective`. This will be stored in `default_test$plot` as a `ggplot2` object that can be manipulated using `ggplot2` functionality. To disable plotting, indicate using the `generate_plot = FALSE` argument.

<ul>

<li>

<h3>Return filtered output</h3>

</li>

</ul>

Alternatively, instead of returning an annotated data frame, `detect_damage` can return the filtered count matrix `filter_threshold = TRUE`. Here, just as above, filtering is done according to a threshold for the estimated level of damage, specified using the `filter_threshold` parameter. By default, `filter_threshold = 0.75`. Essentially, this provides a means of automating the filtering process shown above.

``` r
# Perform damage detection & filtering according to non-default threshold
filtered_test <- detect_damage(
  count_matrix = test_counts, 
  filter_counts = TRUE,
  filter_threshold = 0.7
)

# View output
head(filtered_test)

dim(filtered_test)
#
```

To explore these parameters and more, you can visit the package vignette
 or via the function help page,

``` r
?detect_damage()
```

## License

`DamageDetective` is made available for public use through the [GNU AGPL-3.0](https://opensource.org/licenses/AGPL-3.0) license.

## Authors

**Alicen Henning**\
Stellenbosch University, Cape Town, South Africa\
Bioinformatics and Computational Biology

This work was done under the supervision of Prof Marlo Möller, Prof Gian van der Spuy, and Prof André Loxton.

## References

- McGinnis, C. S., Murrow, L. M., & Gartner, Z. J. (2019). DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. *Cell Systems, 8*(4), 329-337.e4. [https://doi.org/10.1016/j.cels.2019.03.003](https://doi.org/10.1016/j.cels.2019.03.003)

------------------------------------------------------------------------
