data("test_counts", package = "DamageDetective")

test_that("Detect incorrect organism", {

  expect_error(get_organism_indices(count_matrix = test_counts,
                                    organism = "Mmus"))

})

test_that("Handle non-default organism", {

  user_organism <- list(mito_pattern = "^MT-",
                     ribo_pattern = "^(RPS|RPL)",
                     nuclear = c("NEAT1","XIST", "MALAT1"))

  pattern_idx <- get_organism_indices(count_matrix = test_counts,
                                      organism = user_organism)

  expect_gt(length(pattern_idx$mito_idx), 0)

})
