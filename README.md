# DamageDetective <img src="man/figures/logo.svg" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cosimameyer/overviewR/actions) ![Build Status](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/build.yml/badge.svg)

<!-- badges: end -->

## Content

[Description](#description) \| [Installation](#installation) \| [Quick start](#quick-start) \| [Authors](#authors) \| [License](#license) \| [References](#references)

## Description

Low-quality cells in single-cell RNA sequencing (scRNA-seq) include those that succumb to stress before being isolated for sequencing. Damaged cells distort analyses by exhibiting altered mRNA profiles compared to their viable counterparts, making their detection and removal from scRNA-seq data a standard practice [1](#references).

Current approaches identify damaged cells based on deviations in quality control metrics influenced by the well-accepted concept of damage: loss of plasma membrane integrity [2](#references). However, this approach assumes that all viable cells follow similar distributions for these metrics—an assumption that does not always prove effective in heterogeneous samples [3](#references), [4](#references). More recent methods improve upon this by analyzing cells at a population level, grouping those with similar distributions before detecting damaged cells within them. However, this assumes that all distinct populations represent true cells, meaning that, if abundant, damaged cells risk misclassification. Ultimately, the filtering decisions of current approaches rely more on statistical definitions of deviation within a dataset than on biological definitions of damage.

`DamageDetective` takes a different approach inspired by a technique popularised by `DoubletFinder` [5](#references) for detecting doublets: instead of detecting damaged cells by comparing them to each other, detect damaged cells by comparing them to artificially damaged versions of themselves. Using the same principle of damage—loss of plasma membrane integrity—it simulates artificial damage by selectively depleting cytoplasmic RNA. By comparing the expression profiles of real cells to artificial cells in reduced-dimensional space, `DamageDetective` estimates the degree of damage in the form of a score from 0 to 1, where 1 represents the highest RNA loss and greatest likelihood of damage. This provides an intuitive scale for filtering that is comparable across cell types, samples, and experiments and is driven directly by biological definitions of damage.

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
# [1] 10009   500
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

To explore these parameters and more, you can visit the package [vignette](./vignettes/DamageDetective.html)
 or read directly within `R` via the function help page,

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

1.  Luecken, Malte D, and Fabian J Theis. 2019. “Current Best Practices in Single-Cell RNA-Seq Analysis: A Tutorial.” *Molecular Systems Biology* 15 (6): e8746. [https://doi.org/10.15252/msb.20188746](https://doi.org/10.15252/msb.20188746)
2.  Amezquita, Robert A., Aaron T. L. Lun, Etienne Becht, Vince J. Carey, Lindsay N. Carpp, Ludwig Geistlinger, Federico Marini, et al. 2020. “Orchestrating Single-Cell Analysis with Bioconductor.” *Nature Methods* 17 (2): 137–45. [https://doi.org/10.1038/s41592-019-0654-x](https://doi.org/10.1038/s41592-019-0654-x)
3.  Osorio, Daniel, and James J Cai. 2021. “Systematic Determination of the Mitochondrial Proportion in Human and Mice Tissues for Single-Cell RNA-Sequencing Data Quality Control.” *Bioinformatics* 37 (7): 963–67. [https://doi.org/10.1093/bioinformatics/btaa751](https://doi.org/10.1093/bioinformatics/btaa751)
4.  Montserrat-Ayuso, Tomàs, and Anna Esteve-Codina. 2024. “Revealing the Prevalence of Suboptimal Cells and Organs in Reference Cell Atlases: An Imperative for Enhanced Quality Control.” *BioRxiv*, April. [https://doi.org/10.1101/2024.04.18.590104](https://doi.org/10.1101/2024.04.18.590104)
5.  McGinnis, C. S., Murrow, L. M., & Gartner, Z. J. (2019). DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. *Cell Systems, 8*(4), 329-337.e4. [https://doi.org/10.1016/j.cels.2019.03.003](https://doi.org/10.1016/j.cels.2019.03.003)

------------------------------------------------------------------------
