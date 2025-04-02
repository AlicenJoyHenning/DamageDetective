# DamageDetective <img src="man/figures/logo.svg" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cosimameyer/overviewR/actions) ![Build Status](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/build.yml/badge.svg)

<!-- badges: end -->

## Content

[Description](#description) \| [Installation](#installation) \| [Quick start](#quick-start) \| [Authors](#authors) \| [License](#license) \| [References](#references)

Jump to [website](https://alicenjoyhenning.github.io/DamageDetective/ "https://alicenjoyhenning.github.io/DamageDetective/")

## Description

Damaged cells are a class of low-quality artifact targeted during the quality control (QC) of single-cell RNA-seq (scRNA-seq) data. These are cells that, for unknown reasons, succumbed to stress before being sequenced and as a result are associated with gene expression data that fails to describe the cells in their true, viable states.

Current approaches detect damage according to deviations in cell-level quality control metrics. These approaches assume all viable cells follow similar distributions across QC metrics, an assumption that rarely holds in heterogeneous samples. More recent methods address this by analysing cells at a population level, isolating cells with similar distributions before assessing finer deviations. However, this assumes that all distinct distributions are associated with true cells, meaning abundant damage is at risk of being undetected. Ultimately, the filtering decisions of current approaches are controlled more by statistical definitions of deviation than biological definitions of damage.

`DamageDetective` takes a different approach inspired by `DoubletFinder`$^1$, a high performing tool for doublet QC, another low quality scRNA-seq artifact. Here, rather than detecting damage by measuring the extent to which cells deviate from each other, damage is detected by measuring the extent to which cells deviate from artificially damaged profiles of themselves. `DamageDetective` estimates the damage severity of true cells as a score from 0 to 1, providing an intuitive, reproducible scale for filtering that is standardised across cell types, samples, and experiments.

<br>

## Installation

`DamageDetective` can be installed from CRAN,

``` r
install.packages('DamageDetective')
```

Or directly from GitHub,

```         
library(devtools)
devtools::install_github("AlicenJoyHenning/DamageDetective", build_vignettes = TRUE)
```

To verify installation, run the following to see if you can view the package vignette and the function help pages,

``` r
library(DamageDetective)
help(package = "DamageDetective")
```

<br> <br>

## Quick start

The demonstrations below can be followed immediately after loading the package and serve as a test to ensure all is running smoothly. For more detailed examples and explanations, please refer to the package vignette.


### Prepare input

Damage detection is carried out by the `detect_damage` function. This requires only a count matrix to run. We will use the dummy count matrix provided by `DamageDetective`, `test_counts`, a subset of the [(kotliarov-pbmc-2020)](10.1038/s41591-020-0769-8%5D) PBMC dataset provided in the `scRNAseq`$^2$ package.

``` r
library(DamageDetective)
data("test_counts", package = "DamageDetective")
dim(test_counts)
# [1] 32738   500
```

<br>

### Select parameters for damage detection

#### `ribosome_penalty`

While `detect_damage` requires only a count matrix as input, there are optional parameters that adjust the implementation of the function. Of these, we recommend `ribosome_penalty` be adjusted for each input dataset. This is done automatically using the `select_penalty` function. This requires the count matrix as input and will output a numeric of the optimal penalty.

``` r
penalty <- select_penalty(
  count_matrix = test_counts,
  seed = 7
)

# View penalty
penalty
# [1] 1e-05
```

#### `filter_threshold`

`DamageDetective` does not provide binary classifications, "cell" or "damage", it provides a score from 0 to 1 indicating the estimated extent of damage in the cell. This is taken directly from the extent of cytoplasmic RNA loss, a proportion ranging from 0 to 1, of the set of artificial cells which a true cell has the greatest proximity to.

In other words, the user must choose a value between 0 and 1 to determine the level of damage, or proportion of RNA loss, above which cells will be excluded. By default, `DamageDetective` offers the threshold of `0.7`. Values greater than `0.7` reflect more permissive filtering while those closer to `0` reflect more stringent filtering. We recommend the default for all cases but suggest that if adjustments are made, they are informed by inspecting the output `detect_damage` plots, `generate_plot = TRUE`.

The remaining parameters for `detect_damage` can be explored in the [function help page](https://alicenjoyhenning.github.io/DamageDetective/docs/reference/detect_damage.html) or in the [introduction article](https://alicenjoyhenning.github.io/DamageDetective/articles/detection-vignette.html) online.

<br>

### Run damage detection

Damage detection can then be run using the count matrix and selected parameters as inputs. Below we have specified `filter_counts` parameter to `TRUE`. This will use the default `filter_threshold` to identify cells for removal and return the filtered count matrix as output.

``` r
# Perform damage detection
detection_results <- detect_damage(
  count_matrix = test_counts,
  ribosome_penalty = penalty,
  filter_counts = TRUE,
  seed = 7
)
# Simulating cells between 1e-05 and 0.08 RNA loss...
# Simulating cells between 0.1 and 0.3 RNA loss...
# Simulating cells between 0.3 and 0.5 RNA loss...
# Simulating cells between 0.5 and 0.7 RNA loss...
# Simulating cells between 0.7 and 0.9 RNA loss...
# Computing pANN...

# View the resulting count matrix
dim(detection_results$output)
# [1] 32738   470
```

> Note, the above assumes the data is of human origin, see `organism` parameter.

::: {style="text-align: center;"}
![Plot showing the output of test_counts with cells coloured according to the estimated level of damage](man/figures/plot.svg)
:::

<br>

Alternatively, if `filter_counts` is set to `FALSE`, a data frame will be given as output containing the damage scores for each barcode. This is provided for interest to the user if they wish to interact with the `DamageDetective` results more directly. From here, a user can filter their data as done by `filter_counts`.

``` r
# Perform damage detection
detection_results <- detect_damage(
  count_matrix = test_counts,
  ribosome_penalty = penalty,
  filter_counts = FALSE,
  seed = 7
)
# Simulating cells between 1e-05 and 0.08 RNA loss...
# Simulating cells between 0.1 and 0.3 RNA loss...
# Simulating cells between 0.3 and 0.5 RNA loss...
# Simulating cells between 0.5 and 0.7 RNA loss...
# Simulating cells between 0.7 and 0.9 RNA loss...
# Computing pANN...

# View output
print(head(detection_results$output), row.names = FALSE)
#                     Cells DamageDetective
# TCTGGAAAGCCCAACC_H1B2ln6      0.03826609
# CCGTTCATCGTGGGAA_H1B2ln2      0.48333333
# CTTCTCTTCAGCCTAA_H1B2ln1      0.05217739
# GGATTACAGGGATGGG_H1B2ln1      0.01044348
# TCTATTGTCTGGTATG_H1B2ln2      0.02435478
# ACGGGTCAGACAAGCC_H1B2ln6      0.28888889

# Filter matrix 
undamaged_cells <- filter(detection_results$output, DamageDetective < 0.7)
filtered_matrix <- test_counts[, undamaged_cells]
dim(filtered_matrix)
# [1] 32738   470
```

## License

`DamageDetective` is made available for public use through the [GNU AGPL-3.0](https://opensource.org/licenses/AGPL-3.0) license.

## Authors

**Alicen Henning**\
Stellenbosch University, Cape Town, South Africa\
Bioinformatics and Computational Biology

This work was done under the supervision of Prof Marlo Möller, Prof Gian van der Spuy, and Prof André Loxton.

## References

1.  McGinnis, C. S., Murrow, L. M., & Gartner, Z. J. (2019). DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. *Cell Systems, 8*(4), 329-337.e4. <https://doi.org/10.1016/j.cels.2019.03.003>

2.  Risso D, Cole M (2024). *scRNAseq: Collection of Public Single-Cell RNA-Seq Datasets*. <doi:10.18129/B9.bioc.scRNAseq> <https://doi.org/10.18129/B9.bioc.scRNAseq>, R package version 2.20.0, <https://bioconductor.org/packages/scRNAseq>.

------------------------------------------------------------------------
