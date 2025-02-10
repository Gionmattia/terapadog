# Unit tests for get_FCs
library(testthat)
library(terapadog)

# Tests for detect_separator()
test_that("detect_separator works correctly for .csv files", {
  expect_equal(terapadog:::detect_separator("example.csv"), ",")
})

test_that("detect_separator works correctly for .tsv files", {
  expect_equal(terapadog:::detect_separator("example.tsv"), "\t")
})

test_that("detect_separator throws an error for unsupported file formats", {
  expect_error(terapadog:::detect_separator("example.txt"),
               "Unsupported file format. Please provide a .csv or .tsv file.")
})

test_that("detect_separator throws an error for files with no extension", {
  expect_error(terapadog:::detect_separator("example"),
               "Unsupported file format. Please provide a .csv or .tsv file.")
})


# Test for empty or NULL data frames

test_that("check_input_df should throw an error if given a NULL input", {
  expect_error(
    terapadog:::check_input_df(NULL),
    "One of the input dataframes is NULL. Check again inputs")
  })

test_that("check_input_df throws an error for an empty dataframe", {
  df1 <- data.frame()
  expect_error(
    terapadog:::check_input_df(df1),
    fixed = TRUE,
    "One of input dataframes is empty (no rows or columns)")
})

# Tests for check_matching_colnames()

test_that("check_matching_colnames works for matching column names", {
  df1 <- data.frame(A = 1, B = 2, C = 3)
  df2 <- data.frame(C = 3, B = 2, A = 1) # Same column names, different order
  expect_no_error(terapadog:::check_matching_colnames(df1, df2))
})

test_that("check_matching_colnames throws an error for mismatched column names", {
  df1 <- data.frame(A = 1, B = 2, C = 3)
  df2 <- data.frame(A = 1, B = 2, D = 3) # D instead of C
  expect_error(
    terapadog:::check_matching_colnames(df1, df2),
    "SampleNames do not match between the RNA and RIBO/POLY count files."
  )
})

test_that("check_matching_colnames throws an error for partially matching column names", {
  df1 <- data.frame(A = 1, B = 2)
  df2 <- data.frame(A = 1, B = 2, C = 3) # Extra column C
  expect_error(
    terapadog:::check_matching_colnames(df1, df2),
    "SampleNames do not match between the RNA and RIBO/POLY count files."
  )
})

test_that("check_matching_colnames works for dataframes with reordered columns", {
  df1 <- data.frame(A = 1, B = 2, C = 3)
  df2 <- data.frame(C = 3, A = 1, B = 2) # Same columns, different order
  expect_no_error(terapadog:::check_matching_colnames(df1, df2))
})

# Tests for check_value_range()

test_that("check_value_range throws an error for values between 0 and 1", {
  df <- data.frame(A = runif(10, 0, 1), B = runif(10, 0, 1))  # Create dataframe with values in [0, 1]
  expect_error(
    terapadog:::check_value_range(df),
    "Value range in the data is between 0 and 1. Data submitted must be RAW COUNT DATA, not normalised or scaled!"
  )
})

