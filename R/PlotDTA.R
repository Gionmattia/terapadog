# R/plotDTA.R
#'
#' This function will plot an interactive html plot of the results of get_FCs.R
#' That is to say, a plot of the genes undergoing translational regulation,
#' coloured by RegMode.
#' Genes whose RegMode was Undeterminable or Undetermined are omitted.
#' @param FC_results A matrix containing the counts from RNA and RIBO
#' samples.
#' @param path A string, pointing to where to save the html plot.
#' @export

plotDTA <- function(FC_results, path) {
  # Filters out omitted RegModes
  df <- FC_results %>%
    filter(!("RegMode" %in% c("Undeterminable", "Undetermined")))

  # Defines custom colour palette for each Regmode
  custom_colors <- c("Buffered" = "#7ACAFF", "Exclusive" = "#FE939F",
                     "Forwarded" = "#73C3B0", "Intensified" = "#FFD07E",
                     "No Change" = "gray")

  # Create the plot

  plot <- plot_ly(
    data = df,
    x = ~RNA_FC,
    y = ~RIBO_FC,
    text = ~paste("GeneID:", Identifier),
    color = ~RegMode,
    colors = custom_colors,
    type = 'scatter',
    mode = 'markers',
    hovertemplate = "%{text}<extra></extra>"
  ) %>%
    layout(
      title = "Differential Translation Plot",
      xaxis = list(title = "RNA Counts Fold Change"),
      yaxis = list(title = "RIBO Counts Fold Change"),
      plot_bgcolor = "#FFFFFF",
      paper_bgcolor = "#FFFFFF"
    )

  saveWidget(plot, path)
  return(plot)

}
