# Retrieve dummy data ----
data("test_counts", package = "DamageDetective")
test_counts <- test_counts[, 1:100]

# Input data of correct form ----
test_that("verify simulate_counts inputs", {
  # Non-sparse matrix
  test_counts_dense <- Matrix::as.matrix(test_counts)
  expect_warning(simulate_counts(
    count_matrix = test_counts_dense,
    damage_proportion = 0.1,
    generate_plot = FALSE,
    seed = 42))

  # No count matrix
  expect_error(simulate_counts(
    damage_proportion = 0.1, generate_plot = FALSE, seed = 42))

  # Incorrect input for count matrix
  expect_error(simulate_counts(
    count_matrix = "string", damage_proportion = 0.1,
    generate_plot = FALSE, seed = 42
  ))

  # Damage proportion out of range
  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = 1.2
  ))

  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = -1.2
  ))

  # beta_shape_parameters invalid
  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.1,
    beta_shape_parameters = 1
  ))

  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.1,
    beta_shape_parameters = c(-1, 4)
  ))

  # Target damage out of range
  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.1,
    target_damage = 0.7
  ))

  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.1,
    target_damage = c(0, 1.2)
  ))

  # Steepness not one of the defaults
  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.1,
    distribution_steepness = "flat"
  ))

  # Distributions not one of the defaults
  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.1,
    damage_distribution = "right"
  ))

  # Ribosomal penalty out of range
  expect_error(simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.1,
    ribosome_penalty = "high"
  ))

})


# Output data of correct form ----
test_that("output of simulate_counts of correct form", {

  simulate_output <- simulate_counts(
      count_matrix = test_counts, damage_proportion = 0.1,
      generate_plot = FALSE, seed = 7)

  # Dimensions unchanged by simulation
  expect_equal(dim(simulate_output$matrix), dim(test_counts))

  # Summary data frame of correct form (not necessarily order)
  expect_true(setequal(
    colnames(simulate_output$qc_summary),
    c('Cell', 'Damaged_Level', 'Original_Features', 'New_Features',
      'Original_MitoProp', 'New_MitoProp', 'Original_RiboProp',
      'New_RiboProp')
  ))

})


# Function logic testing ----

test_that("Simulating no perturbation", {

  simulate_output <- simulate_counts(
    count_matrix = test_counts, damage_proportion = 0,
    generate_plot = FALSE, seed = 7)

  expect_identical(simulate_output$qc_summary$Original_MitoProp,
               simulate_output$qc_summary$New_MitoProp)

  expect_identical(simulate_output$qc_summary$Original_RiboProp,
               simulate_output$qc_summary$New_RiboProp)

})

test_that("Simulating complete perturbation", {

  simulate_output <- simulate_counts(
    count_matrix = test_counts,
    damage_proportion = 1,
    target_damage = c(0.5, 1),
    generate_plot = FALSE, seed = 7)

  expect_false(any(simulate_output$qc_summary$Original_MitoProp ==
                     simulate_output$qc_summary$New_MitoProp))

  expect_false(any(simulate_output$qc_summary$Original_RiboProp ==
               simulate_output$qc_summary$New_RiboProp))

})


# Reproducibility testing ----

test_that("Testing simulate_counts seed", {

  counts_1 <- simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.05,
    generate_plot = FALSE, seed = 777)

  counts_2 <- simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.05,
    generate_plot = FALSE, seed = 777)

  counts_3 <- simulate_counts(
    count_matrix = test_counts, damage_proportion = 0.05,
    generate_plot = FALSE, seed = 1)

  # If same seed, output should be identical
  expect_identical(counts_1$matrix, counts_2$matrix)
  expect_identical(counts_1$qc_summary, counts_2$qc_summary)

  # If different seed, output should not be identical
  expect_false(all(counts_1$matrix == counts_3$matrix))
  expect_false(all(counts_1$qc_summary == counts_3$qc_summary))

})



