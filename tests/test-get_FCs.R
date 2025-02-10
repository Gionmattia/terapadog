# Unit tests for get_FCs
library(testthat)
library(terapadog)

test_that("get_FCs raises an error for NULL inputs", {
  expect_error(get_FCs(NULL, data.frame()), "expression.data or exp_de cannot be NULL.")
  expect_error(get_FCs(matrix(), NULL), "expression.data or exp_de cannot be NULL.")
})

# Input Check - exp_de required columns
test_that("Input validation for exp_de", {
  expression.data <- matrix(1:9, nrow = 3, dimnames = list(NULL, c("S1", "S2", "S3")))
  # Check if required columns are present
  expect_error(get_FCs(expression.data, data.frame(Sample = c("S1", "S2", "S3"))), "exp_de is missing required columns: .*")
})

test_that("get_FCs raises an error if exp_de is empty after Samples with Group = NA are removed", {
            # create a dataframe with Group = NA for everything in exp_de.
            exp_de <- data.frame(
              SampleID = c("S1", "S2"),
              Group = c(NA, NA),
              SeqType = c("RNA", "RIBO")
            )
            # create a df where SampleID is the column name
            expression.data <- data.frame(
              S1 = c(1,2),
              S2 = c(2,1)
            )

  expect_error(get_FCs(expression.data, exp_de), "exp_de cannot be empty.")
})

test_that("get_FCs raises an error if samples in expression.data do not match with the ones from exp_de", {
  # Create a matrix with 3 columns
  expression.data <- matrix(1:9, nrow = 3, dimnames = list(NULL, c("S1", "S2", "S3")))

  # Create a data frame with 2 rows instead of 3
  exp_de <- data.frame(
    SampleID = c("S1", "S2"),
    Group = c("A", "B"),
    SeqType = c("RNA", "RIBO")
  )

  # Expect an error
  expect_error(
    get_FCs(expression.data, exp_de),
    "The SampleID (colnames(expression.data)) must match the SampleID values in exp_de.",
    fixed = TRUE
  )
})


test_that("get_FCs checks for 2 comaprison groups (2 levels in Group) and multiple
          blocks (more than 1 level for Block) in exp_de after NA removal", {
  expression.data <- matrix(1:9, nrow = 3, dimnames = list(NULL, c("S1", "S2", "S3")))

  # Case 1: Single comparison group in exp_de
  exp_de_single_group <- data.frame(
    SampleID = c("S1", "S2", "S3"),
    Group = c("A", "A", "A"),
    SeqType = c("RNA", "RIBO", "RNA"),
    Block = c(1, 1, 1)
  )
  expect_error(
    get_FCs(expression.data, exp_de_single_group),
    "exp_de must contain two levels for Group - The differential analysis requires to compare two groups only, not more nor less."
  )

  # Case 2: Single block with paired = TRUE
  exp_de_single_block <- data.frame(
    SampleID = c("S1", "S2", "S3"),
    Group = c("A", "B", "A"),
    SeqType = c("RNA", "RIBO", "RNA"),
    Block = c(1, 1, 1)
  )
  expect_error(
    get_FCs(expression.data, exp_de_single_block, paired = TRUE),
    "When paired = TRUE, exp_de must contain more than one level for Block."
  )
})

