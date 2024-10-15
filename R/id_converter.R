# R/id_converter.R

#' Convert the human gene identifier (hgnc_symbol or ensembl_gene_id) to
#' entrezgene_id format for the analysis.
#' @aliases id_converter
#' @importFrom biomaRt useMart getBM
#' @importFrom stats setNames
#' @importFrom utils write.table
#' @param esetm A matrix with the gene count values and whose rownames are the
#' gene Ids (gene symbol or ensembl gene ID).
#' @param id_type A string representing the type of ID given as input. Must be
#' either hgnc_symbol or ensembl_gene_id.
#' @return A matrix with gene IDs in the entrezgene_id format. Also provides a
#' report on the duplicated mappings (conversion_report.txt) in the working dir.
#' @examples
#' # To showcase thisl internal function, a small example is made.
#' gene_ids <- c('ENSG00000103197', 'ENSG00000008710', 'ENSG00000167964'
#' , 'ENSG00000167964')
#' esetm <- matrix(c(
#' 2.5, 3.1, 5.2, 0.1,
#' 4.1, 2.9, 6.3, 0.5,
#' 1.5, 3.7, 4.8, 0.1), nrow = 4, byrow = FALSE)
#' rownames(esetm) <- gene_ids
#' colnames(esetm) <- c("Sample 1", "Sample 2", "Sample 3")
#' # Call the function
#' esetm <- id_converter(esetm, "ensembl_gene_id")
#' print(head(esetm))
#' @export
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

  # Extract the duplicated mappings (one entrezgene_id could map to more ids)
  duplicates_mask <- duplicated(gene_mapping$entrezgene_id) | duplicated(gene_mapping$entrezgene_id, fromLast = TRUE)
  duplicated_rows <- gene_mapping[duplicates_mask, ]
  duplicated_rows <- duplicated_rows[!is.na(duplicated_rows$entrezgene_id), ]
  # Write them to a report for transparency
  report_message <- paste0(
    "These are the input Ids mapping to the same entrezgene_id. ",
    "NA values (no mapping) have been removed. Be aware that the counts from ",
    "multiple mappings have been aggregated (summed) together under the same ",
    "entrezgene_id for the purpose of the 'terapadog' analysis. ",
    nrow(duplicated_rows), " duplicates found:"
  )

  # Write the explanation to a file
  report_path <- paste0(getwd(),"/conversion_report.txt")
  writeLines(report_message, report_path)
  # Append the table to the file
  write.table(duplicated_rows, report_path, append = TRUE,
              sep = ",", row.names = FALSE, col.names = TRUE)

  # Create a named vector for mapping
  id_to_entrez <- stats::setNames(gene_mapping$entrezgene_id, gene_mapping[[id_type]])

  # Map the IDs in esetm to Entrez IDs
  entrez_ids <- id_to_entrez[rownames(esetm)]

  # Update the row names of esetm
  rownames(esetm) <- entrez_ids

  # Remove rows with NA Entrez IDs
  esetm <- esetm[!is.na(rownames(esetm)), ]

  # Aggregate the multiple mappings (duplicates) under the same entrezgene_id
  esetm <- rowsum(esetm, group = rownames(esetm))

  # Remove duplicates (if any) maybe not necessary
  # esetm <- esetm[!duplicated(rownames(esetm)), ]

  return(esetm)
}

