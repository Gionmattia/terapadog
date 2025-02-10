# Unit tests for plotDTA
library(testthat)
library(terapadog)

# Input Check - FC_results
test_that("Input validation for FC_results", {
  # Check if NULL
  expect_error(plotDTA(NULL), "FC_results must be a dataframe")
  # Check if required columns are present
  expect_error(plotDTA(data.frame()), "FC_results is missing required columns: .*")
})

# Checks output
test_that("plotDTA filters out 'Undeterminable' and 'Undetermined' RegModes before plotting", {
  df <- data.frame(
    Identifier = c("Gene A", "Gene B", "Gene C"),
    RegMode = c("Buffered", "Undeterminable", "No Change"),
    RNA_FC = c(-0.4, 0.3, 0.1),
    RIBO_FC = c(0.2, -0.1, 0.3)
  )
  result <- plotDTA(df)

  # Extract the data used for the plot
  plot_data <- result$x$data[[1]]$text

  # Ensure "Undeterminable" and "Undetermined" are not in the plot
  expect_false(any(grepl("Undeterminable", plot_data)))
  expect_false(any(grepl("Undetermined", plot_data)))

  # tests on saving conditions
  # Test 1: Error is raised when directory does not exist
  expect_error(
    plotDTA(df, save_plot = TRUE, path = "fake/path"),
    regexp = "The specified directory does not exist: .*"
  )

  # Test 2: No error is raised when save_plot = FALSE
  expect_silent(
    plotDTA(df, save_plot = FALSE, path = "fake/path")
  )

  # Test 3: No error is raised when the directory exists (in this case, temp_dir)
  expect_silent(
    plotDTA(df, save_plot = TRUE)
  )
})

