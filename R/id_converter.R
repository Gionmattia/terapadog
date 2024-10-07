# Conversion to Entrez gene ID
# This could be used as a preprocessng helper function...
# or maybe other to the Entrez... this has some potential!

# R/id_converter.R

#' Convert the human gene identifier (hgnc_symbol or ensembl_gene_id) to
#' entrezgene_id format for the analysis.
#' @importFrom biomaRt useMart getBM
#' @importFrom stats setNames
#' @param esetm A dataframe with the gene count values.
#' @param id_type A string representing the type of ID given as input. Must be
#' either hgnc_symbol or ensembl_gene_id.
#' @return A dataframe with gene IDs in the entrezgene_id format.
#' @examples
#' # Internal function, the following code cannot be run from here.
#' \dontrun{
#' # To showcase this internal function, only one file (rna_counts) is used.
#' path <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
#' esetm <- read.table(path, header = TRUE, sep = "\t")
#' # Call the function
#' converted_esetm <- id_converter(esetm, "ensembl_gene_id")
#' print(head(esetm))
#' }
#' @keywords internal
#'

id_converter <- function(esetm, id_type) {

  if (!(id_type %in% c("hgnc_symbol", "ensembl_gene_id"))) {
    stop("Error: the conversion script only works with
         'hgnc_symbol' or 'ensembl_gene_id'")
  }

  # Use the Ensembl database
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  # The IDs to map are the row names of esetm
  ids <- rownames(esetm)

  # Get the mapping between the given IDs and Entrez Gene IDs
  gene_mapping <- biomaRt::getBM(
    filters = id_type,
    attributes = c(id_type, "entrezgene_id"),
    values = ids,
    mart = mart
  )

  # Create a named vector for mapping
  id_to_entrez <- stats::setNames(gene_mapping$entrezgene_id, gene_mapping[[id_type]])

  # Map the IDs in esetm to Entrez IDs
  entrez_ids <- id_to_entrez[rownames(esetm)]

  # Update the row names of esetm
  rownames(esetm) <- entrez_ids

  # Remove rows with NA Entrez IDs
  esetm <- esetm[!is.na(rownames(esetm)), ]

  # Remove duplicates (if any)
  esetm <- esetm[!duplicated(rownames(esetm)), ]

  return(esetm)
}

