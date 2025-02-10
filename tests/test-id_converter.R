# Unit tests for get_FCs
library(testthat)
library(terapadog)

test_that("unsupported geneID types are not processed", {
  gene_ids <- c('ENSG00000103197', 'ENSG00000008710', 'ENSG00000167964'
                , 'ENSG00000167964')
                esetm <- matrix(c(
                2.5, 3.1, 5.2, 0.1,
                4.1, 2.9, 6.3, 0.5,
                1.5, 3.7, 4.8, 0.1), nrow = 4, byrow = FALSE)
                rownames(esetm) <- gene_ids
                colnames(esetm) <- c("Sample 1", "Sample 2", "Sample 3")

  expect_error(id_converter(esetm, "a_wrong_id"), "Error: the conversion script only works with
         'hgnc_symbol' or 'ensembl_gene_id'")
})

# Checks the existance of terapadog's gene mapping file. This is a .csv containing a
# conversion table between geneIDs, to use if the Ensembl server is down. More info in the extdata_info.md file
test_that("Terapadog's gene mapping file (to use if ensembl is down) exists within the package", {
  data_path <- system.file("extdata", "hsapiens_GRCh38.p14_gene_mapping.csv.gz", package = "terapadog")
  expect_true(file.exists(data_path), info = "Local gene mapping file is missing. ")
})
