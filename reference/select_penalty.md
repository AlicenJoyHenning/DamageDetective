# select_penalty

Recommended prerequisite function to detect_damage() that estimates the
ideal `ribosome_penalty` value for the input data.

## Usage

``` r
select_penalty(
  count_matrix,
  organism = "Hsap",
  mito_quantile = 0.75,
  penalty_range = c(0.1, 1),
  penalty_step = 0.05,
  max_penalty_trials = 14,
  target_damage = c(0.2, 0.99),
  damage_distribution = "right_skewed",
  distribution_steepness = "steep",
  beta_shape_parameters = NULL,
  stability_limit = 3,
  damage_proportion = 0.15,
  annotated_celltypes = FALSE,
  return_output = "penalty",
  ribosome_penalty = NULL,
  seed = 7,
  verbose = TRUE
)
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

- mito_quantile:

  Numeric specifying below what proportion of mitochondrial content
  cells are used for sampling for simulation.

  - Default is 0.75, meaning only cells with less than 0.75 proportion
    of mitochondrial counts are sampled for simulated.

- penalty_range:

  Numerical vector of length 2 specifying the lower and upper limit of
  values tested for the ribosomal penalty.

  - Default is c(0.1, 1).

- penalty_step:

  Numeric specifying the value added to each increment of penalty
  tested.

  - Default is 0.05.

- max_penalty_trials:

  Numeric specifying the maximum number of iterations for the ribosomal
  penalty value.

  - Default is 14.

- target_damage:

  Numeric vector specifying the upper and lower range of the level of
  damage that will be introduced.

  Here, damage refers to the amount of cytoplasmic RNA lost by a cell
  where values closer to 1 indicate more loss and therefore more heavily
  damaged cells.

  - Default is c(0.65, 1)

- damage_distribution:

  String specifying whether the distribution of damage levels among the
  damaged cells should be shifted towards the upper or lower range of
  damage specified in 'target_damage' or follow a symmetric distribution
  between them. There are three valid options:

  - "right_skewed"

  - "left_skewed"

  - "symmetric"

  - Default is "left_skewed"

- distribution_steepness:

  String specifying how concentrated the spread of damaged cells are
  about the mean of the target distribution specified in
  'target_damage'. Here, an increase in steepness manifests in a more
  apparent skewness.There are three valid options:

  - "shallow"

  - "moderate"

  - "steep"

  - Default is "steep"

- beta_shape_parameters:

  Numeric vector that allows for the shape parameters of the beta
  distribution to defined explicitly. This offers greater flexibility
  than allowed by the 'damage_distribution' and 'distribution_steepness'
  parameters and will override the defaults they offer.

  - Default is 'NULL'

- stability_limit:

  Numeric specifying the number of additional iterations allotted after
  the median minimum distance of the artificial cells to the true cells
  is greater than the previous minimum distance.

  The idea here is that if a higher penalty is not causing an
  improvement in the output, there is little need to continue testing
  with larger penalties.

  - Default is 3.

- damage_proportion:

  Numeric describing what proportion of the input data should be altered
  to resemble damaged data.

  - Must range between 0 and 1.

- annotated_celltypes:

  Boolean specifying whether input matrix has cell type information
  stored.

  - Default is FALSE

- return_output:

  String specifying what form the output of the function should take
  where the options are either,

  - "penalty"

  - "full"

  "Penalty" will return only the ribosomal penalty that resulted in the
  best performance (the smallest median distance between artificial and
  true cells). While "full" will return the ideal ribosomal penalty and
  the median distance between artificial and true cells for each penalty
  tested. This allows insight into how the penalty was selected.

  - Default is "penalty".

- ribosome_penalty:

  Numeric specifying the factor by which the probability of loosing a
  transcript from a ribosomal gene is multiplied by. Here, values closer
  to 0 represent a greater penalty.

  - Default is 1.

- seed:

  Numeric specifying the random seed to ensure reproducibility of the
  function's output. Setting a seed ensures that the random sampling and
  perturbation processes produce the same results when the function is
  run multiple times with the same input data and parameters.

  - Default is 7.

- verbose:

  Boolean specifying whether messages and function progress should be
  displayed in the console.

  - Default is TRUE.

## Value

Numeric representing the ideal ribosomal penalty for an input dataset.

## Details

Based on observations of true single cell data, we find that ribosomal
RNA loss occurs less frequently than expected based on abundance alone.
To adjust for this, the probability scores of ribosomal gene loss are
multiplied by a numerical value (`ribosome_penalty`) between 0 and 1.
Lower values (closer to zero) better approximate true data, with a
default of 0.01, though this can often be greatly refined for the input
data.

Refinement follows a similar workflow to detect_damage(), but rather
than evaluating the similarity of true cells to sets of artificial cells
to infer their level of damage, we evaluate the similarity of artificial
cells to true cells to infer the effectiveness of their approximation to
true data. This is calculated using the distance to the nearest true
cell (dTNN) taken for each artificial cell found using the Euclidean
distance matrix. The median dTNN is computed iteratively until
stabilization or a worsening trend. The ideal `ribosomal_penalty` is
then selected as that which generated the lowest dTNN.

## Examples

``` r
data("test_counts", package = "DamageDetective")

penalty <- select_penalty(
 count_matrix = test_counts,
 stability_limit = 1,
 max_penalty_trials = 1,
 seed = 7
)
#> Testing penalty of 0.1...
#> Maximum penalty trials reached (1). Stopping.
```
