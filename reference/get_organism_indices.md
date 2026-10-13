# Retrieve genes corresponding to the organism of interest

Retrieve genes corresponding to the organism of interest

## Usage

``` r
get_organism_indices(count_matrix, organism)
```

## Arguments

- count_matrix:

  Matrix or dgCMatrix containing the counts from single cell RNA
  sequencing data.

- organism:

  String specifying the organism of origin of the input data where there
  are two standard options,

  - "Hsap"

  - "Mmus"

  If a user wishes to use a non-standard organism they must input a list
  containing strings for the patterns to match mitochondrial and
  ribosomal genes of the organism. If available, nuclear-encoded genes
  that are likely retained in the nucleus, such as in nuclear speckles,
  must also be specified. An example for humans is below,

  - organism = list(mito_pattern = "^MT-", ribo_pattern = "^(RPS\|RPL)",
    nuclear = c("NEAT1","XIST", "MALAT1"))

  - Default is "Hsap"

## Value

List containing the indices of the count matrix corresponding to
mitochondrial, non-mitochondrial, and ribosomal gene sets.
