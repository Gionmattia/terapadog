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

# Debug: testing get_FCs

resultsss <- get_FCs(expression.data, exp_de)

# Debug: testing terapadog itself

#load gene.indices
load("~/Desktop/terapadog/data/gene_identifier_set.Rdata")

is_paired <- FALSE

res <- terapadog(esetm = expression.data, exp_de = exp_de, paired = is_paired,
                 NI = 1000, Nmin = 0)

test <- as.character(unlist(gene.indices))
