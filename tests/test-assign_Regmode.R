library(testthat)
library(terapadog)

# Test how function handles empty dataframe
test_that("assign_Regmode handles empty dataframe", {
  # Create an empty dataframe with the required columns
  empty_df <- data.frame(
    Identifier = character(0),
    log2FoldChange = numeric(0),
    padj = numeric(0),
    RIBO_FC = numeric(0),
    RIBO_padj = numeric(0),
    RNA_FC = numeric(0),
    RNA_padj = numeric(0)
  )
  # test a dataframe with no genes (nrows = 0)
  expect_error(terapadog:::assign_Regmode(empty_df), "Input dataframe has no rows.", fixed = TRUE)
})

# Test behavior if columns are missing
test_that("assign_Regmode handles missing dataframe columns", {
  expect_error(terapadog:::assign_Regmode(data.frame()), "Input is missing required columns: .*")
  })

# Test for "Undeterminable" RegMode -> where at least one padj is NA
test_that("assign_Regmode assigns 'Undeterminable' correctly", {
  df <-  data.frame(
    Identifier = c("NA_test1", "NA_test2", "NA_test3", "NA_test4"),
    log2FoldChange = c(0,0,0,0),
    padj = c(NA, 0.02, 0.03, NA),
    RIBO_FC = c(0,0,0, 0),
    RIBO_padj = c(0.3, NA, 0.03, NA),
    RNA_FC = c(0,0,0,0),
    RNA_padj = c(0.3, 0.2, NA, NA)
  )
  result <- terapadog:::assign_Regmode(df)
  expect_true(all(result$RegMode == "Undeterminable"))
})

# Test for "Undetermined" RegMode -> as per defined by Chotani et al 2019
test_that("assign_Regmode assigns 'Undetermined' correctly", {
  df <-  data.frame(
    Identifier = c("Un_test1", "Un_test2", "Un_test3"),
    log2FoldChange = c(0,0,0),
    padj = c(0.6, 0.6, 0.01),
    RIBO_FC = c(0,0,0),
    RIBO_padj = c(0.6, 0.01, 0.06),
    RNA_FC = c(0,0,0),
    RNA_padj = c(0.01, 0.06, 0.06)
  )
  result <- terapadog:::assign_Regmode(df)
  expect_true(all(result$RegMode == "Undetermined"))
})

# Test for the "No Change" Regmode
test_that("assign_Regmode assigns 'No Change' correctly", {
  df <-  data.frame(
    Identifier = c("NC_test"),
    log2FoldChange = c(0),
    padj = c(0.6),
    RIBO_FC = c(0),
    RIBO_padj = c(0.6),
    RNA_FC = c(0),
    RNA_padj = c(0.06)
  )
  result <- terapadog:::assign_Regmode(df)
  expect_true(all(result$RegMode == "No Change"))
})

# Test for the "Intensified" Regmode. This condition requires all padj to be <0.05
# and log2FoldChange * RNA_FC > 1  (aka same directionality of change)
test_that("assign_Regmode assigns the 'Intensified' cases correctly", {
  df <-  data.frame(
    Identifier = c("Intensified_test1", "Intensified_test2"),
    log2FoldChange = c(1, -1),
    padj = c(0.01, 0.01),
    RIBO_FC = c(0, 0),
    RIBO_padj = c(0.01, 0.01),
    RNA_FC = c(1, -1),
    RNA_padj = c(0.01, 0.01)
  )
  result <- terapadog:::assign_Regmode(df)
  expect_true(all(result$RegMode == "Intensified"))
})

# Test for the "Buffered" Regmode.
test_that("assign_Regmode assigns the 'Buffered' cases correctly", {
  df <-  data.frame(
    Identifier = c("Buffered_test1", "Buffered_test2", "Buffered_test3"),
    log2FoldChange = c(1, -1, 1),
    padj = c(0.01, 0.01, 0.01),
    RIBO_FC = c(1, 1, 1),
    RIBO_padj = c(0.01, 0.01, 0.6),
    RNA_FC = c(-1, 1, 1),
    RNA_padj = c(0.01, 0.01, 0.01)
  )
  result <- terapadog:::assign_Regmode(df)
  expect_true(all(result$RegMode == "Buffered"))
})

# Test for the "Exclusive" Regmode.
test_that("assign_Regmode assigns the 'Exclusive' cases correctly", {
  df <-  data.frame(
    Identifier = c("Exclusive_test"),
    log2FoldChange = c(1),
    padj = c(0.01),
    RIBO_FC = c(1),
    RIBO_padj = c(0.01),
    RNA_FC = c(0),
    RNA_padj = c(0.6)
  )
  result <- terapadog:::assign_Regmode(df)
  expect_true(all(result$RegMode == "Exclusive"))
})

# Test for the "Forwarded" Regmode.
test_that("assign_Regmode assigns the 'Forwarded' cases correctly", {
  df <-  data.frame(
    Identifier = c("Forwarded_test"),
    log2FoldChange = c(0),
    padj = c(0.6),
    RIBO_FC = c(1),
    RIBO_padj = c(0.01),
    RNA_FC = c(1),
    RNA_padj = c(0.01)
  )
  result <- terapadog:::assign_Regmode(df)
  expect_true(all(result$RegMode == "Forwarded"))
})
