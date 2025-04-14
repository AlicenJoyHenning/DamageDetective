# Retrieve dummy data ----
data("test_counts", package = "DamageDetective")
test_counts <- test_counts[, 1:100]

simulate_output <- simulate_counts(count_matrix = test_counts,
                                   damage_proportion = 0.1, seed = 777)

detection_output <- data.frame(
  features = runif(100, 0, 2000),
  mt.prop = runif(100, 0, 1),
  rb.prop = runif(100, 0, 1),
  DamageDetective = runif(100, 0, 1)
)

# Outputs of correct form ----
test_that("ggplot2 output for plot_altered_counts", {
  plot <- plot_altered_counts(simulate_output$qc_summary)
  expect_s3_class(plot, "ggplot")
})

test_that("ggplot2 output for plot_unaltered_counts", {
  plot <- plot_unaltered_counts(simulate_output$qc_summary)
  expect_s3_class(plot, "ggplot")
})

test_that("ggplot2 output for plot_simulation_output", {
  plot <- plot_simulation_outcome(simulate_output$qc_summary)
  expect_s3_class(plot, "ggplot")
})

test_that("ggplot2 output for plot_ribosomal_opnalty", {
  plot <- plot_ribosomal_penalty(simulate_output$qc_summary)
  expect_s3_class(plot, "ggplot")
})

test_that("ggplot2 output for plot_detection_output", {
  plot <- plot_detection_outcome(detection_output)
  expect_s3_class(plot, "ggplot")
})
