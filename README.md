# DamageDetective

[![R-CMD-check](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/R-CMD-check.yaml) ![Build Status](https://github.com/AlicenJoyHenning/DamageDetective/actions/workflows/build.yml/badge.svg)


R package for understanding and managing the presence of damaged cells in single cell RNA sequencing data.

## Documentation Content

[Description](#description) \| [Installation](#installation) \| [Quickstart](#quickstart) \| [Authors](#authors) \| [License](#license) \| [References](#references)

## Description

Single-cell RNA sequencing (scRNA-seq) is a well-established technique in the era of next-generation sequencing. However, it relies heavily on the quality of upstream pre-processing, an extensive collection of steps to go from raw sequencing files to a cell-type annotated count matrix. Despite being a fundamental step of the pre-processing workflow, cell-level quality control, particularly the removal of damaged cells, is neglected.

We propose that the difficulty in detecting damaged cells stems from the difficulty in defining them. Unlike other low-quality cell artifacts, damage is not a binary classification but exists on a nonlinear, stochastic spectrum where not all will qualify for filtering. Without a sure definition of damage, it is difficult to know how to best go about detecting them.

`DamageDetective` provides a computational solution by exploiting a traditional characteristic of damaged cells: the loss of plasma membrane integrity. In single cell data, this manifests as the increased likelihood of RNA leakage in damaged cells and that can be simulated through the probabilistic loss of cytoplasmic RNA. Using a simulation framework, `DamageDetective` predicts which cells of an input matrix are most likely damaged by comparing them to their simulated damaged counterparts. 



## Installation

`DamageDetective` can be installed from CRAN using,

``` r
install.packages('DamageDetective')
```

Alternatively, if you want to use the latest development version, you can install the pacakge directly from GitHub using,

``` r
library(devtools)
devtools::install_github("AlicenJoyHenning/DamageDetective",ref='devel')
```

::: {style="gray"}
To verify installation, run the following to see if you can view the package vignette and the function help page,

``` r
library(DamageDetective)
help(package = "DamageDetective")
?DamageDetective()
```
:::

------------------------------------------------------------------------

## Quick start

The demonstrations below can be followed immediately after loading the package and serve as a "test run" to ensure all is running smoothly. For more advanced function descriptions and usage examples please refer to the [package vignettes](link).

Begin by loading the dummy data provided by the package, `test_counts`, a PBMC dataset containing 500 cells and 10009 genes.

``` r
data("test_counts", package = "DamageDetective")
dim(test_counts)
# [1] 10009   500
```

### Damaged cell simulation

#### \- Default

A core task of `DamageDetective` is to predict how cells from existing scRNA-seq datasets might appear if they had experienced a certain degree of damage. Here, damage is modeled by the loss of cytoplasmic RNA where cells with great RNA loss are assumed to be extensively damaged, while those with minimal loss are considered largely intact.

This achieved using the `simulate_counts` function. While `simulate_counts` is used internally by the `detect_damage` function, it has been made available to the public as an exploratory tool providing a flexible framework to understand damage as it manifests in single cell data. It may also prove useful in the generation of annotated ground truth datasets.

Using an input `count_matrix` and a `target_proportion` of damage, `simulate_counts` will return a count matrix with a randomly selected target proportion of damaged cells with altered gene expression profiles. Below, we are introducing damage to 25 % of cells from the input data, 125 of the `test_counts`.

The altered counts are stored in the `matrix` slot of the output.

``` r
# Run the simulation
simulated_counts <- simulate_counts(
  count_matrix = test_counts, 
  damage_proportion = 0.25
)

# View the output 
dim(simulated_counts) 
# [1] 10009   500 # counts remain unchanged
head(simulated_counts)
```

In addition to the altered count matrix is a data frame containing quality control statistics for the cells before and after simulation stored in the `qc_summary` slot of the output. This can be helpful for understanding how the damaged cells are changing, particularly according to the level of damage it was assigned.

``` r
head(simulated_counts$qc_summary) 
```

By default, the output is visualised in a plot grid returned by `simulate_counts`. The plots show the distribution of cells according to quality control metrics before alteration, displayed in the top row, and after alteration, in the bottom row. This is stored in the `plot` slot of the output as a `ggplot2` object. To disable plotting, indicate using the `generate_plot = FALSE` argument.

#### \- Flexibility

But the power of the `simulate_counts` function comes in the flexibility offered by its input parameters mainly,

-   `target_damage` specifying the upper and lower range of damage that will be introduced across the selected cells.
-   `damage_distribution` specifying the distribution of damage introduced across the cells within the `target_damage` range.

For example, to generate damaged cells that are only extensively damaged we can specify a narrow target, between 70 and 99 % of RNA loss, that is shifted heavily towards the upper limit, 99 % loss, we can say,

``` r
# Run the simulation
simulated_counts <- simulate_counts(
  count_matrix = test_counts, 
  damage_proportion = 0.25, 
  target_damage = c(0.7, 0.99),
  damage_distribution = "right_skewed"
)
```

Alternatively, to generate a wide range of damaged cells that are mostly very lightly damaged we can specify a wider target, between 0.01 and 99 % of RNA loss, that is shifted towards the lower limit, 0.01 % loss, we can say,

``` r
# Run the simulation
simulated_counts <- simulate_counts(
  count_matrix = test_counts, 
  damage_proportion = 0.25, 
  target_damage = c(0.01, 0.99),
  damage_distribution = "left_skewed"
)
```

To explore these parameters and more, you can visit the package [vignette](link) or read directly within `R` via the function help page,

``` r
?simulate_counts()
```

### Damaged cell detection

#### \- Default

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

#### \- Return filtered output

Alternatively, instead of returning an annotated data frame, `detect_damage` can return the filtered count matrix `filter_threshold = TRUE`. Here, just as above, filtering is done according to a threshold for the estimated level of damage, specified using the `filter_threshold` parameter. By default, `filter_threshold = 0.75`.

Essentially, this provides a means of automating the filtering process shown above.

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

To explore these parameters and more, you can visit the package [vignette](link) or read directly within `R` via the function help page,

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

------------------------------------------------------------------------
