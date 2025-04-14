# Retrieve dummy data ----
data("test_counts", package = "DamageDetective")
test_counts <- test_counts[, 1:100]

# Input data of correct form ----
test_that("verify select_penalty inputs", {
  # Inputs out of acceptable range
  expect_error(select_penalty(count_matrix = test_counts, mito_quantile = 12))
  expect_error(select_penalty(count_matrix = test_counts, penalty_step = 2))
  expect_error(select_penalty(count_matrix = test_counts,
                              penalty_range = c(0.1, 1.2)))
  expect_error(select_penalty(count_matrix = test_counts,
                              max_penalty_trials = 5000))
  expect_error(select_penalty(count_matrix = test_counts,
                              stability_limit = 5000))


  # Inputs of incorrect format
  expect_error(select_penalty(count_matrix = test_counts, mito_quantile = "1"))
  expect_error(select_penalty(count_matrix = test_counts, penalty_range = 1))
  expect_error(select_penalty(count_matrix = test_counts, penalty_step = "1"))
  expect_error(select_penalty(count_matrix = test_counts,
                              return_output = TRUE))

})

# Output data of correct form ----
test_that("select_penalty output correctness", {

  penalty_output <- select_penalty(count_matrix = test_counts,
                                   max_penalty_trials = 2,
                                   seed = 7, return_output = "full",
                                   verbose = FALSE)

  expect_type(select_penalty(count_matrix = test_counts,
                             max_penalty_trials = 2,
                            seed = 7, verbose = FALSE), "double")


  expect_true(all(c('Penalty', 'Global_mean') %in%
            colnames(penalty_output$penalty_results))
    )

})

# Logic of function ----
test_that("select_penalty logic", {

  penalty_output <- select_penalty(count_matrix = test_counts,
                                   max_penalty_trials = 2,
                                   seed = 7, return_output = "full",
                                   verbose = FALSE)

  # Within range
  expect_gte(penalty_output$selected_penalty, 0)
  expect_lte(penalty_output$selected_penalty, 1)

  # Correct selection
  expect_true(min(penalty_output$penalty_results) ==
              penalty_output$selected_penalty)

})

# Reproducibility testing ----
test_that("Testing select_penalty seed", {

  penalty_1 <- select_penalty(
    count_matrix = test_counts, max_penalty_trials = 4,
    verbose = FALSE, return_output = "full", seed = 777)

  penalty_2 <- select_penalty(
    count_matrix = test_counts, max_penalty_trials = 4,
    verbose = FALSE, return_output = "full", seed = 777)

  penalty_3 <- select_penalty(
    count_matrix = test_counts, max_penalty_trials = 4,
    verbose = FALSE, return_output = "full", seed = 778)

  expect_identical(penalty_1$penalty_results, penalty_2$penalty_results)
  expect_false(isTRUE(all.equal(penalty_1$penalty_results,
                                penalty_3$penalty_results)))

})
