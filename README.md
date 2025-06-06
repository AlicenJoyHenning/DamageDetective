# DamageDetective <img src="man/figures/logo.svg" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cosimameyer/overviewR/actions) ![Build Status](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/build.yml/badge.svg) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![CRAN status](https://www.r-pkg.org/badges/version/DamageDetective)](https://CRAN.R-project.org/package=DamageDetective)

<!-- badges: end -->

## Content

[Description](#description) \| [Installation](#installation) \| [Quick start](#quick-start) \| [Contribute](#contribute) \| [Authors](#authors) \| [License](#license) \| [References](#references)

## Description

[**Jump to the DamageDetective website**](https://alicenjoyhenning.github.io/DamageDetective/)

Damaged cells are an artifact of single-cell RNA sequencing (scRNA-seq) formed when cells succumb to stress before being sequenced. As a result, the gene expression data captured does not reflect biologically viable cells and introduces technical variability that is indistinguishable from functionally relevant variability. Filtering these cells is a standard task in scRNA-seq quality control (QC), though lacks standardisation in practice.

The majority of approaches filter damaged cells according to deviations in cell-level QC metrics. This outlier-based detection implicitly assumes viable cells follow unimodal distributions across QC metrics, where deviation is synonymous with damage. This assumption falters in the context of heterogeneous data and risks introducing filtering bias related to cell type abundance. Recent methods address this by defining damage within distinct distributions, representing cell populations, independently. This, however, assumes all distinct distributions are associated with viable cell populations and risks leaving abundant damage undetected and ultimately misclassified. 

`DamageDetective` takes a different approach, rather than detecting damage by measuring the extent to which cells deviate from one another, it measures the extent to which cells deviate from artificially damaged profiles of themselves, created through simulating cytoplasmic RNA escape–a characteristic of damage resulting from the loss of plasma membrane integrity. This is inspired by the approach of `DoubletFinder`—a high-performing tool of another prominent scRNA-seq artifact. 

Like `DoubletFinder`, `DamageDetective` uses principal component analysis to compute the proximity of true cells to artificial cells. This is calculated as a proportion (pANN) of a cell's nearest neighbours that are of artificial origin, reflecting the likelihood that the cell has experienced the same cytoplasmic RNA loss as its artificial neighbours, i.e., is damaged. This score, ranging from 0 to 1, provides an intuitive scale for filtering that is standardised across cell types, sample origin, and experimental design.

<br>

## Installation

Install `DamageDetective` from CRAN (R \>= 4.4.0),

``` r
install.packages('DamageDetective')
```

Or the latest development version on GitHub (R \>= 3.5.0),

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

This demonstration can be followed immediately after loading the package using the internal dummy dataset. For examples with true data and more detailed explanations, please refer to the package articles [website](https://alicenjoyhenning.github.io/DamageDetective/). 

### Prepare input

Damage detection is carried out by the `detect_damage` function that accepts count matrices, `Seurat` or `SingleCellExperiment` objects, or alignment files ([package tutorials](https://alicenjoyhenning.github.io/DamageDetective/articles/detection-vignette.html)) as input. We will demonstrate using a dummy count matrix, `test_counts`, a subset of the [(kotliarov-pbmc-2020)](10.1038/s41591-020-0769-8%5D) PBMC dataset provided in the `scRNAseq` package.

``` r
library(DamageDetective)
library(Matrix)
data("test_counts", package = "DamageDetective")
dim(test_counts)
```
> Expected outcome, 
> ```R 
> [1] 32738   500
> ```


<br>

### Select parameters for damage detection

#### `ribosome_penalty`

While `detect_damage` requires only a count matrix as input, additional parameters control aspects of the function's computations. Of these, we recommend `ribosome_penalty` be adjusted for each dataset using the `select_penalty` function as shown below,

``` r
penalty <- select_penalty(count_matrix = test_counts)
penalty
```
> Expected outcome, 
> ```R 
> Testing penalty of 0.1...
> Testing penalty of 0.15...
> Testing penalty of 0.2...
> Testing penalty of 0.25...
> Stopping early: dTNN is no longer improving.
>
> 0.1
> ```


#### `filter_threshold`

`DamageDetective` performs filtering using the proximity scores according to a threshold. By default, `DamageDetective` offers the threshold of `0.5` where values greater than `0.5` reflect more permissive filtering and values closer to `0` reflect more stringent filtering. We recommend the default, but suggest that if adjustments are made, they are informed by the output `detect_damage` plots, `generate_plot = TRUE`.

<br>

### Run damage detection

Damage detection is run as shown below, using the count matrix and ribosomal penalty as inputs. Below, we have additionally specified for `filter_counts` parameter to be `TRUE`. This will use the default `filter_threshold` to detect damaged cells for removal and return the filtered count matrix that can be used immediately afterwards for the remainder of pre-processing. Though implemented in R, `DamageDetective` provides output that is platform-agnostic and can be integrated into any existing single-cell analysis workflow.

``` r
# Perform damage detection
detection_results <- detect_damage(
  count_matrix = test_counts,
  ribosome_penalty = penalty,
  filter_counts = TRUE
)

# View the resulting count matrix
dim(detection_results$output)
```
> Expected outcome, 
> ```R 
> Clustering cells...
> Simulating damage...
> Computing pANN...
>
>  32738   461
> ```

<br>

Alternatively, if `filter_counts` is set to `FALSE`, a data frame will be given as output containing the damage scores for each cell. This is provided for the user if they wish to interact with the `DamageDetective` results directly. From here, a user can filter their data manually, as is done by `filter_counts=TRUE` automatically.

``` r
# Perform damage detection
detection_results <- detect_damage(
  count_matrix = test_counts,
  ribosome_penalty = penalty,
  filter_counts = FALSE,
  seed = 7
)

# View output
print(head(detection_results$output), row.names = FALSE)

# Filter matrix 
undamaged_cells <- subset(detection_results$output, DamageDetective < 0.7)
filtered_matrix <- test_counts[, undamaged_cells$Cells]
dim(filtered_matrix)
```
> Expected outcome, 
> ```R 
> Clustering cells...
> Simulating damage...
> Computing pANN...
>
>                      Cells DamageDetective
> TCTGGAAAGCCCAACC_H1B2ln6               0
> CCGTTCATCGTGGGAA_H1B2ln2               0
> CTTCTCTTCAGCCTAA_H1B2ln1               0
> GGATTACAGGGATGGG_H1B2ln1               0
> TCTATTGTCTGGTATG_H1B2ln2               0
> ACGGGTCAGACAAGCC_H1B2ln6               0
> 
> 32738   461
> ```

## Contribute

We are committed to the improvement of `DamageDetective` and encourage users to report any bugs or difficulties they encounter. Contributions that refine or challenge the assumptions and heuristics used to detect damaged cells are also welcome. Please reach out via the maintainer's email listed in the `DESCRIPTION` file or start a public discussion [![Issue](https://img.shields.io/badge/Issues-blue?style=flat&logo=github)](https://github.com/AlicenJoyHenning/DamageDetective/issues). 

## License 

`DamageDetective` is made available for public use through the [GNU AGPL-3.0](https://opensource.org/license/agpl-v3)

## Authors 

**Alicen Henning**\
Stellenbosch University, Cape Town, South Africa\
Bioinformatics and Computational Biology

This work was done under the supervision of Prof Marlo Möller, Prof Gian van der Spuy, and Prof André Loxton.

## References 

- McGinnis, C. S., Murrow, L. M., & Gartner, Z. J. (2019). DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. *Cell Systems, 8*(4), 329-337.e4. <https://doi.org/10.1016/j.cels.2019.03.003>

- Risso D, Cole M (2024). *scRNAseq: Collection of Public Single-Cell RNA-Seq Datasets*. <doi:10.18129/B9.bioc.scRNAseq> <https://doi.org/10.18129/B9.bioc.scRNAseq>, R package version 2.20.0, <https://bioconductor.org/packages/scRNAseq>.
