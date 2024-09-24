# R/get_FCs.R

#'
#' This function execute the Differential Translation Analysis on its own
#' using DeltaTE.
#' The output is a dataframe with the FC in mRNA counts, RIBO counts or TE
#' between the conditions in exam.
#' @param expression.data A matrix containing the counts from RNA and RIBO
#' samples.
#' @param exp_de A dataframe containing information regarding the samples.
#' It has number of rows equal to the columns of esetm.
#' @param paired Logical. Default is false. Set to TRUE if the experiment
#' has paired samples in its design.
#' @export
#'
#'

get_FCs <- function(expression.data, exp_de, paired = FALSE) {

  # DeltaTE cannot manage NA among the levels of a factor, so I need to remove them! (code same as above)
  na_samples <- exp_de$SampleID[is.na(exp_de$Group)]
  # Remove the NA samples from exp_de
  exp_de <- exp_de[!exp_de$SampleID %in% na_samples, ]
  # Remove NA samples from expression.data
  expression.data <- expression.data[, !colnames(expression.data) %in% na_samples, drop = FALSE]

  if (paired == TRUE) {
    design_TE <- ~ Block + Group + SeqType + Group:SeqType
    design_R <- ~ Block + Group
  }
  else {
    design_TE <- ~ Group + SeqType+ Group:SeqType
    design_R <- ~ Group
  }

  # run DeltaTE to get the TE linear model
  ddsMat <- DESeqDataSetFromMatrix(
    countData = expression.data, # Ribo_counts and rna_counts should have been provided as a single dataframe already.
    colData = exp_de, # This is where we give DeltaTe the info
    design = design_TE
  )

  # Run DESeq
  ddsMat <- suppressMessages(DESeq(ddsMat))

  # Extract the fold changes
  res <- results(ddsMat, name = "Groupd.SeqTypeRIBO")

  # FC analysis for RNA counts:
  # Create filter using exp_de to find samples with SeqType "RNA"
  rna_samples <- exp_de %>%
    filter("SeqType" == "RNA") %>%
    select("SampleID")

# Extract Sample IDs as a vector
  rna_sample_ids <- rna_samples$SampleID

# Select columns from expression.data
  expression.data.rna <- expression.data[, rna_sample_ids]

# Execute the analysis for RNA counts
  ddsMat_rna <- DESeqDataSetFromMatrix(
    countData = expression.data.rna,
    colData=exp_de[which(exp_de$SeqType == "RNA"),],
    design = design_R)

  ddsMat_rna <- DESeq(ddsMat_rna)
  res_rna <- results(ddsMat_rna, name="Group_d_vs_c") # Check name of result! we encode it as "d" and "c"
  res_rna <- lfcShrink(ddsMat_rna,coef="Group_d_vs_c", res=res_rna)

  # FC analsysis for RIBO counts:
  # Create filter using exp_de to find samples with SeqType  "RIBO"
  ribo_data <- exp_de %>%
    filter("SeqType" == "RIBO") %>%
    select("SampleID")

  # Extract Sample IDs as a vector
  ribo_sample_ids <- ribo_data$SampleID

  # Select columns from expression.data
  expression.data.ribo <- expression.data[, ribo_sample_ids]

  # Getting FC for Ribo
  ddsMat_ribo <- DESeqDataSetFromMatrix(
    countData = expression.data.ribo,
    colData = exp_de[which(exp_de$SeqType == "RIBO"),],
    design = design_R)

  ddsMat_ribo <- DESeq(ddsMat_ribo)
  res_ribo <- results(ddsMat_ribo,name="Group_d_vs_c")
  res_ribo <- lfcShrink(ddsMat_ribo,coef="Group_d_vs_c", res=res_ribo)

  # convert results to dataframe
  res$Identifier <- rownames(res)
  res_combined <- as.data.frame(res)

  res_rna_slim <- as.data.frame(res_rna) %>%
    select(RNA_FC = "log2FoldChange", RNA_padj = "padj")

  res_ribo_slim <- as.data.frame(res_ribo) %>%
    select(RIBO_FC = "log2FoldChange", RIBO_padj = "padj")

  res_rna_slim$Identifier <- rownames(res_rna_slim)
  res_ribo_slim$Identifier <- rownames(res_ribo_slim)

  res_combined <- res_combined %>%
    left_join(res_ribo_slim, by = "Identifier") %>%
    left_join(res_rna_slim, by = "Identifier")

  res_combined <- assign_Regmode(res_combined)

  # Reorder columns
  col_order <- c("Identifier", colnames(res_combined)[1:ncol(res_combined)])

return(res_combined[, col_order])
}
