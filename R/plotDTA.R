# R/plotDTA.R

#' This function will plot an interactive html plot of the results of get_FCs.R
#' That is to say, a plot of the genes undergoing translational regulation,
#' coloured by RegMode.
#' Genes whose RegMode was Undeterminable or Undetermined are omitted.
#' @importFrom dplyr %>% filter sym
#' @importFrom plotly plot_ly layout
#' @importFrom htmlwidgets saveWidget
#' @param FC_results A matrix containing the counts from RNA and RIBO
#' samples.
#' @param path A string, pointing to where to save the html plot. If none is
#' given, the plot will be saved to a temporary directory.
#' @return An interactive html plot.
#' @examples
#' # Creates a mock dataframe for this demonstration
#' df <- data.frame(
#'   Identifier = c("Gene A", "Gene B", "Gene C", "Gene D"),
#'   RegMode = c("Buffered", "Exclusive", "Undeterminable", "No Change"),
#'   RNA_FC = c(-0.40, -0.5, NA, 0.01),
#'   RIBO_FC = c(0.19, -0.3, 0.8, -0.02)
#' )
#' result <- plotDTA(df)
#' @export
plotDTA <- function(FC_results, path = paste0(tempdir(), "/plot.html")) {
  # Filters out omitted RegModes
  df <- FC_results %>%
    dplyr::filter(!(!!dplyr::sym("RegMode") %in% c("Undeterminable", "Undetermined")))

  # Defines custom colour palette for each Regmode
  custom_colors <- c("Buffered" = "#7ACAFF", "Exclusive" = "#FE939F",
                     "Forwarded" = "#73C3B0", "Intensified" = "#FFD07E",
                     "No Change" = "gray")

  # Create the plot

  plot <- plotly::plot_ly(
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

  htmlwidgets::saveWidget(plot, path)
  return(plot)

}
