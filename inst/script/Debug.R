# Debug: everything to get the data ready for terapadog
path_to_RNA_counts <- "~/Desktop/terapadog/inst/extdata/rna_counts.tsv"
path_to_RIBO_counts <- "~/Desktop/terapadog/inst/extdata/ribo_counts.tsv"
path_to_metadata <- "~/Desktop/terapadog/inst/extdata/sample_info.tsv"
analysis.group.1 <- "1"
analysis.group.2 <- "2"


dats <- prepareTerapadogData(path_to_RNA_counts, path_to_RIBO_counts,
                             path_to_metadata,
                             analysis.group.1, analysis.group.2)

expression.data <- dats$expression.data
exp_de <- dats$exp_de

# Debug: testing terapadog itself

#load gene.indices
#load("~/Desktop/terapadog/data/gene_identifier_set.Rdata")

expression.data <- id_converter(expression.data, "ensembl_gene_id")

is_paired <- FALSE

res <- terapadog2(esetm = expression.data, exp_de = exp_de, paired = is_paired,NI = 1000, Nmin = 0)

# Troubleshooting old code.
res <- keggLink("pathway", "hsa")
a <- data.frame(path = gsub(paste0("path:", "hsa"),"", res),
                gns = gsub(paste0("hsa", ":"),"", names(res)))
gslist <- tapply(a$gns, a$path, function(x) {
  as.character(x)
})
gs.names <- keggList("pathway", "hsa")[paste0("hsa", names(gslist))]
names(gs.names) <- names(gslist)

# So this should solve the issue and allow the use of... "no gene.indices" in the


## OLD VERSION

res <- keggLink("pathway", organism)
a <- data.frame(path = gsub(paste0("path:", organism),
                            "", res), gns = gsub(paste0(organism, ":"),
                                                 "", names(res)))
gslist <- tapply(a$gns, a$path, function(x) {
  as.character(x)
})
gs.names <- keggList("pathway", organism)[paste0("path:",
                                                 organism, names(gslist))]
names(gs.names) <- names(gslist)



#' # Data is available in the "inst/extdata" folder of this package
prepared_data <- prepareTerapadogData("./inst/extdata/rna_counts.tsv", "./inst/extdata/ribo_counts.tsv", "./inst/extdata/sample_info.tsv", "1", "2")
#' # Unpacks the expression.data and exp_de from the output
expression.data <- prepared_data$expression.data
exp_de <- prepared_data$exp_de
#' # For sake of brevity, only the data frame's head will be printed out
print(head(expression.data))
print(head(exp_de))



# Trying get_FCs.R

results <- get_FCs(expression.data, exp_de)
