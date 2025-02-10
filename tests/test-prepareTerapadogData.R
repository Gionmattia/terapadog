library(testthat)
library(terapadog)
# tests/testthat/test-prepareTerapadogData.R

test_that("Fails if RNA counts file does not exist", {
  rna_file <- system.file("extdata", "not_rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv", package = "terapadog")
  meta_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")

  expect_error(
    terapadog:::prepareTerapadogData(
      path_to_RNA_counts = rna_file,
      path_to_RIBO_counts = ribo_file,
      path_to_metadata = meta_file,
      analysis.group.1 = "1",
      analysis.group.2 = "2"
    ),
    regexp = "does not exist"
  )
})

test_that("Fails if RIBO counts file does not exist", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "not_ribo_counts.tsv", package = "terapadog")
  meta_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")

  expect_error(
    terapadog:::prepareTerapadogData(
      path_to_RNA_counts = rna_file,
      path_to_RIBO_counts = ribo_file,
      path_to_metadata = meta_file,
      analysis.group.1 = "1",
      analysis.group.2 = "2"
    ),
    regexp = "does not exist"
  )
})

test_that("Fails if sample info file does not exist", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv", package = "terapadog")
  meta_file <- system.file("extdata", "not_sample_info.tsv", package = "terapadog")

  expect_error(
    terapadog:::prepareTerapadogData(
      path_to_RNA_counts = rna_file,
      path_to_RIBO_counts = ribo_file,
      path_to_metadata = meta_file,
      analysis.group.1 = "1",
      analysis.group.2 = "2"
    ),
    regexp = "does not exist"
  )
})


test_that("Fails if sample_info file is missing required columns", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv", package = "terapadog")
  meta_file <- system.file("extdata", "test_sample_info_wrong.tsv", package = "terapadog")

  expect_error(
    terapadog:::prepareTerapadogData(
      path_to_RNA_counts = rna_file,
      path_to_RIBO_counts = ribo_file,
      path_to_metadata = meta_file,
      analysis.group.1 = "1",
      analysis.group.2 = "2"
    ),
    regexp = "Metadata file is missing required columns: "
  )
})


test_that("Fails if analysis.group.1 and analysis.group.2 are not in metadata Conditions", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv", package = "terapadog")
  meta_file <- system.file("extdata", "test_sample_info_wrong_2.tsv", package = "terapadog")

  expect_error(
    terapadog:::prepareTerapadogData(
      path_to_RNA_counts = rna_file,
      path_to_RIBO_counts = ribo_file,
      path_to_metadata = meta_file,
      analysis.group.1 = "1",
      analysis.group.2 = "2"
    ),
    regexp = "do not match Condition values"
  )
})


test_that("Successful run returns a list with expression.data and exp_de", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv", package = "terapadog")
  meta_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")

  out <- terapadog:::prepareTerapadogData(
    path_to_RNA_counts = rna_file,
    path_to_RIBO_counts = ribo_file,
    path_to_metadata = meta_file,
    analysis.group.1 = "1",
    analysis.group.2 = "2"
  )

  # Check it returns a list with two elements named "expression.data" and "exp_de"
  expect_type(out, "list")
  expect_named(out, c("expression.data", "exp_de"))

  # Optionally check dimensions (should not be necessary, but better safe than sorry)
  expect_equal(ncol(out$expression.data), 16) # Check no samples were left behind
  expect_equal(nrow(out$exp_de), 16) # Checks exp_de has all samples
  expect_true(all(c("SampleID", "SampleName", "SeqType", "Condition", "Group") %in% names(out$exp_de)))
})

