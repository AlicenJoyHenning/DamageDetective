
# DamageDetective <img src="man/figures/logo.svg" align="right" height="139"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/DamageDetective)](https://CRAN.R-project.org/package=DamageDetective)
[![R-CMD-check](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cosimameyer/overviewR/actions)
![Build
Status](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/build.yml/badge.svg)

<!-- badges: end -->

## Content

[Description](#description) \| [Installation](#installation) \| [Quick
start](#quick-start) \| [Authors](#authors) \| [License](#license) \|
[References](#references)

## Description

[**View on DamageDetective
website**](https://alicenjoyhenning.github.io/DamageDetective/)

Damaged cells are an artifact of single-cell RNA sequencing (scRNA-seq)
data formed when cells succumb to stress before, or in the process of,
being sequenced. As a result, the gene expression data captured for
damaged cells fails to reflect viable cell states. Continuing an
analysis with damage compromises its reliability as it is not possible
to differentiate true from technical variability. Filtering damaged
cells is therefore an essential step in scRNA-seq quality control (QC).

Current approaches detect damage according to deviations in cell-level
QC metrics. This outlier-based detection assumes viable cells follow a
unimodal distribution where being deviant is equivalent to being damage.
However, this assumption does not hold in heterogeneous scRNA-seq data. Often, this introduces a
damage filtering bias relative to cell type abundance. More recent approaches address this by detecting damage at a local
scale, isolating cells with similar distributions before addressing
deviations within them. However, this assumes all distinct distributions
are associated with viable cell populations, meaning abundant damage is
at risk of misclassification. Ultimately, the filtering decisions of
current approaches are controlled by statistical definitions of
deviation that do not always align with biological definitions of
damage.

`DamageDetective` takes a different approach, rather than detecting
damage by measuring the extent to which cells deviate from one other, it
measures the extent to which cells deviate from artificially damaged
profiles of themselves. This is inspired by the approach of
`DoubletFinder`, a high performing QC tool for filtering doublets. Using
principal component analysis, the proximity of true cells to sets of
artificial cells with known levels of damage is computed. The damage
level of the set to which a true cell shows the highest proximity is
assigned to the true cell, outputted as a score ranging from 0 to 1.
This provides an intuitive scale for filtering that is standardised
across cell types, sample origin, and experimental design.


> `DamageDetective` runs using platform-agnostic data types in order to integrate seamlessly into any single-cell analysis workflow. However, the popularised R data types for scRNA-seq, `Seurat` and `SingleCellExperiment` objects, as well as alignment output, can also be used ([package
tutorials](https://alicenjoyhenning.github.io/DamageDetective/articles/detection-vignette.html)).

<br>

## Installation

Install `DamageDetective` from CRAN,

``` r
install.packages('DamageDetective')
```

Or through GitHub,

```         
library(devtools)
devtools::install_github("AlicenJoyHenning/DamageDetective", build_vignettes = TRUE)
```

To verify installation, run the following to see if you can view the
package vignette and the function help pages,

``` r
library(DamageDetective)
help(package = "DamageDetective")
```

<br> <br>

## Quick start

The demonstrations below can be followed immediately after loading the
package and serve as a test to ensure all is running smoothly. For more
detailed examples and explanations, please refer to the package articles
available on our
[website](https://alicenjoyhenning.github.io/DamageDetective/).

### Prepare input

Damage detection is carried out by the `detect_damage` function that
requires a count matrix to run. We will use the dummy count matrix
provided by `DamageDetective`, `test_counts`, a subset of the
[(kotliarov-pbmc-2020)](10.1038/s41591-020-0769-8%5D) PBMC dataset
provided in the `scRNAseq`$^2$ package.

``` r
library(DamageDetective)
data("test_counts", package = "DamageDetective")
dim(test_counts)
# [1] 32738   500
```

<br>

### Select parameters for damage detection

#### `ribosome_penalty`

While `detect_damage` requires only a count matrix as input, there are
optional parameters that control the implementation of the function. Of
these, we recommend `ribosome_penalty` be adjusted for each dataset.
This can be done automatically using `select_penalty`, a function that
also requires only the count matrix as input and will output a numeric
of the optimal penalty.

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

`DamageDetective` does not provide binary classifications, "cell" or
"damage", it provides a score from 0 to 1 indicating the estimated
extent of damage in the cell. This is taken directly from the extent of
cytoplasmic RNA loss, a proportion ranging from 0 to 1, of the set of
artificial cells which a true cell has the greatest proximity to.

In other words, the user must choose a value between 0 and 1 to
determine the level of damage, or proportion of RNA loss, above which
cells will be excluded. By default, `DamageDetective` offers the
threshold of `0.7`. Values greater than `0.7` reflect more permissive
filtering while those closer to `0` reflect more stringent filtering. We
recommend the default for all cases but suggest that if adjustments are
made, they are informed by inspecting the output `detect_damage` plots,
`generate_plot = TRUE`.

<br>

### Run damage detection

Damage detection can be run using the count matrix and selected
parameters as inputs. Below we have specified `filter_counts` parameter
to `TRUE`. This will use the default `filter_threshold` to identify
cells for removal and return the filtered count matrix as output.

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

> Note, the above assumes the data is of human origin, see `organism`
> parameter.


![Output of test_counts with cells coloured according to the estimated
level of damage](man/figures/plot.svg)


<br>

Alternatively, if `filter_counts` is set to `FALSE`, a data frame will
be given as output containing the damage scores for each barcode. This
is provided for the user if they wish to interact with the
`DamageDetective` results more directly. From here, a user can filter
their manually data as done by `filter_counts` automatically.

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

`DamageDetective` is made available for public use through the [GNU
AGPL-3.0](https://opensource.org/license/agpl-v3)

## Authors

**Alicen Henning**\
Stellenbosch University, Cape Town, South Africa\
Bioinformatics and Computational Biology

This work was done under the supervision of Prof Marlo Möller, Prof Gian
van der Spuy, and Prof André Loxton.

## References

1.  McGinnis, C. S., Murrow, L. M., & Gartner, Z. J. (2019).
    DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data
    Using Artificial Nearest Neighbors. *Cell Systems, 8*(4),
    329-337.e4. <https://doi.org/10.1016/j.cels.2019.03.003>

2.  Risso D, Cole M (2024). *scRNAseq: Collection of Public Single-Cell
    RNA-Seq Datasets*. <doi:10.18129/B9.bioc.scRNAseq>
    <https://doi.org/10.18129/B9.bioc.scRNAseq>, R package version
    2.20.0, <https://bioconductor.org/packages/scRNAseq>.
