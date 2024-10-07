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
#load("~/Desktop/terapadog/data/gene_identifier_set.Rdata"

# Important!! Need to convert the ID if it's not done before!!
expression.data <- id_converter(expression.data, "ensembl_gene_id")
res <- terapadog(esetm = expression.data, exp_de = exp_de)



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


# creating mockDTA results

mockDTA <- read.table("~/Desktop/terapadog/inst/extdata/mockfile_assignRegmode.tsv", header = TRUE)
mockDTA2 <- assign_Regmode(mockDTA)
write.csv(mockDTA2, "~/Desktop/terapadog/inst/extdata/mock_DTA_results.tsv")

plotDTA(mockDTA2, "~/Desktop/plot_test.html")

count <- sum(mockDTA2$RegMode %in% c("Undeterminable", "Undetermined"))

count

# Downsizing the testing files to reduce the running time of tutorials.

path_to_RNA_counts <- "~/Desktop/terapadog/inst/extdata/rna_counts.tsv"
path_to_RIBO_counts <- "~/Desktop/terapadog/inst/extdata/ribo_counts.tsv"

rna <- read.table(path_to_RNA_counts, sep = "\t")
ribo <- read.table(path_to_RIBO_counts, sep = "\t")

rows_to_remove <- rownames(rna)[1:10000]

# 2. Remove these rows from both df1 and df2
rna <- rna[!rownames(rna) %in% rows_to_remove, ]
ribo <- ribo[!rownames(ribo) %in% rows_to_remove, ]

write.table(rna, "~/Desktop/reduced_rna.tsv", sep = "\t")
write.table(ribo, "~/Desktop/reduced_ribo.tsv", sep = "\t")


# Testing the shortned files

path_to_RNA_counts <- "~/Desktop/reduced_rna.tsv"
path_to_RIBO_counts <- "~/Desktop/reduced_ribo.tsv"

write.table(results, "~/Desktop/reduced_results.tsv", sep = "\t")
plotDTA(results, "~/Desktop/Delete.html")

terapadog(esetm = expression.data, exp_de = exp_de)
