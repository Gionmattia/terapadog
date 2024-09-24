# The following script illustrates how to generate the gene_identifier_set.Rdata
# file within the "data" folder.
# This code has been reported for transparency purposes, it is not intended to
# be run by the user when testing TERAPADOG.

# About libraries: while KEGGREST is part of the TERAPADOG package, biomaRt is
# not. If you wish to run this example, you will have to install it manually.
library(biomaRt)
library(KEGGREST)


# Get the identified genes from the data
# In this example, we are getting the gene IDs from one of the testing datasets
rna_counts <- read.csv("./inst/extdata/rna_counts.tsv", header = TRUE,
                       row.names = 1, sep = "\t")
ensembl_ids <- rownames(rna_counts)


# maps the ENSEMBL gene ids to their ncbi counterpart and removes NAs
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
conversion <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = ensembl)
conversion <- conversion[!is.na(conversion$entrezgene_id), ]


# Maps the Entrez gene IDs to kegg gene IDs
# Obtain the conversion vector from ncbi_geneid to KEGG_id
kegg_ids <- keggConv("hsa", "ncbi-geneid", querySize = 100)
# Merge the converter to kegg_ids. We want a df linking KEGG_Ids to ensembl.
stripped_names <- cleaned_names <- gsub("ncbi-geneid:", "", names(kegg_ids))
kegg_ids_only <- unname(kegg_ids)
first_df <- data.frame(stripped_names, kegg_ids_only)
colnames(first_df)[colnames(first_df) == "stripped_names"] <- "entrezgene_id"


# merge all 3 dfs and remove NA (ensembl GENE ids not mapping to KEGG)
# to get a Final Mapping dataframe.
final_df <- merge(first_df, conversion, by = "entrezgene_id", all.x = TRUE)
final_df <- final_df[!is.na(final_df$ensembl_gene_id), ]


# Ok, now we have a dataframe with a conversion for all the genes that were in
# the data. Time to get the pathways with keggLink
my_Kegg_IDs <- final_df$kegg_ids_only
# Initialize an empty list
kegg_pathways <- list()


# Define a batch size. The process is really long and will be killed many times.
# The following code is made so to allow the user to save their results and
# resubmit the conversion request.
batch_size <- 200

# Loop through the KEGG IDs in batches
for (i in seq(1, length(my_Kegg_IDs), by = batch_size)) {

  # Determine the range for the current batch
  batch_ids <- my_Kegg_IDs[i:min(i + batch_size - 1, length(my_Kegg_IDs))]

  # Process the current batch of KEGG IDs
  for (kegg_id in batch_ids) {

    # Retrieve the pathways associated with the current KEGG ID
    pathways <- keggLink("pathway", kegg_id)

    # Store the retrieved pathways in the list, with the KEGG ID as the name
    kegg_pathways[[kegg_id]] <- pathways

    # Print the current KEGG ID for monitoring progress
    print(kegg_id)
  }

  # Save the progress after each batch
  save(kegg_pathways, file = paste0("kegg_pathways_backup.RData"))

  # Print separator for monitoring progress
  print("Batch complete. Saving progress and sleeping.")
  print("_________")

  # Random sleep between 1 and 5 seconds after processing each batch
  Sys.sleep(sample(1:5, 1))
}

# If the process is interrupted, run the following (and keep doing so until
# all ids have been processed)

load("kegg_pathways_backup.RData")
last_processed_id <- tail(names(kegg_pathways), 1)
start_position <- match(last_processed_id, my_Kegg_IDs) + 1
batch_size <- 100
for (i in seq(start_position, length(my_Kegg_IDs), by = batch_size)) {
  batch_ids <- my_Kegg_IDs[i:min(i + batch_size - 1, length(my_Kegg_IDs))]

  for (kegg_id in batch_ids) {
    pathways <- keggLink("pathway", kegg_id)
    kegg_pathways[[kegg_id]] <- pathways
    print(kegg_id)
  }

  save(kegg_pathways, file = paste0("kegg_pathways_backup.RData"))
  print("Batch complete. Saving progress and sleeping.")
  Sys.sleep(sample(05:30, 1))
}


# Assuming kegg_pathways is your list with KEGG IDs as names and pathways as values
# Initialize an empty list to store the data for the dataframe
pathway_list <- list()
# Loop through the kegg_pathways list to prepare the data
for (kegg_id in names(kegg_pathways)) {
  # Retrieve the pathways associated with the current KEGG ID
  pathways <- kegg_pathways_5[[kegg_id]]

  # Create entries for each pathway
  for (pathway in pathways) {
    if (!is.null(pathway)) {
      pathway_list[[length(pathway_list) + 1]] <- data.frame(
        Source_ID = kegg_id,
        Target_ID = pathway,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine all entries into a single dataframe
pathway_df <- do.call(rbind, pathway_list)

library(dplyr)

# Now you can proceed with further processing, e.g., grouping by pathway
grouped_pathway_df <- pathway_df %>%
  group_by(Target_ID) %>%
  summarise(Genes = paste(Source_ID, collapse = ", "))

# Creates a named list where each entry has the name of a pathway and values are gene IDs.
pathway_list <- setNames(strsplit(grouped_pathway_df$Genes, ", "),
                         grouped_pathway_df$Target_ID)

# Convert the gene IDs to ENSEMBL gene ID.
# Initialize an empty list to store the results
ensembl_list <- list()

# Iterate through the named list
for (pathway in names(pathway_list)) {
  # Get the vector of gene IDs for this pathway
  gene_ids <- pathway_list[[pathway]]
  print(gene_ids)
  # Retrieve the corresponding ENSEMBL IDs from the conversion dataframe
  ensembl_ids <- final_df$ensembl_gene_id[match(gene_ids, final_df$kegg_ids_only)]

  # Store the ENSEMBL IDs in the result list with the pathway as the name
  ensembl_list[[pathway]] <- ensembl_ids
}




## Additional bit: converting gene.indices to the gene_set_list
gene_identifier_set <- lapply(gene.indices, function(gene_ids) {
  rownames(rna_counts)[gene_ids]
})

save(ensembl_list, file = paste0("gene_identifier_set.RData"))
