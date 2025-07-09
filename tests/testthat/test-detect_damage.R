# Retrieve dummy data ----
data("test_counts", package = "DamageDetective")
test_counts <- test_counts[, 1:100]

# Input data of correct form ----
test_that("verify detect_damage inputs", {

  expect_error(detect_damage(
    count_matrix = test_counts, filter_threshold = 1.2, seed = 777
  ))

  expect_error(detect_damage(
    count_matrix = 1.2, seed = 777
  ))


})

# Output data of correct form ----
test_that("correctness of detect_damage filter output", {

  filter_output <- detect_damage(
    count_matrix = test_counts,
    filter_counts = TRUE,
    generate_plot = FALSE,
    seed = 777, verbose = FALSE)

  # Did filtering occur?
  expect_lte(dim(filter_output)[2], dim(test_counts)[2])

})

test_that("correctness of detect_damage unfiltered output", {

  unfiltered_output <- detect_damage(
    count_matrix = test_counts,
    generate_plot = FALSE,
    seed = 777, verbose = FALSE)

  # Are all cells retained?
  expect_equal(dim(unfiltered_output)[1], dim(test_counts)[2])

  # Are the correct columns present?
  expect_true(all(c('DamageDetective') %in%
                    colnames(unfiltered_output))
  )

  # Are damage scores within valid range?
  expect_true(all(unfiltered_output$DamageDetective <= 1))
  expect_true(all(unfiltered_output$DamageDetective >= 0))

})

# Reproducibility testing ----

test_that("Testing detect_damage seed", {
  scores_1 <- detect_damage(
    count_matrix = test_counts,
    verbose = FALSE, generate_plot = FALSE,
    filter_counts = FALSE, seed = 777)

  scores_2 <- detect_damage(
    count_matrix = test_counts,
    verbose = FALSE, generate_plot = FALSE,
    filter_counts = FALSE, seed = 777)

  scores_3 <- detect_damage(
    count_matrix = test_counts,
    verbose = FALSE, generate_plot = FALSE,
    filter_counts = FALSE, seed = 778)

  expect_identical(scores_1$DamageDetective, scores_2$DamageDetective)
  expect_false(isTRUE(all.equal(scores_1$DamageDetective,
                                scores_3$DamageDetective)))
})
