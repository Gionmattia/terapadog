# R/preprocessing_helpers.R

#' Detect the separator used in a file
#'
#' This function determines whether the file uses a comma or tab as a separator
#' based on the file extension (.csv or .tsv).
#'
#' @param path A string representing the file path.
#' @return A string representing the separator (comma or tab).
#' @keywords internal
#'
detect_separator <- function (path) {
  if (grepl("\\.csv$", path)) {
    separator <- ","
  } else if (grepl("\\.tsv$", path)) {
    separator <- "\t"
  } else {
    stop("Unsupported file format. Please provide a .csv or .tsv file.")
  }
  return(separator)
}



#' Check if two data frames have the same column names
#'
#' This function verifies that two data frames (RIBO and RNA) contain
#' the same set of column names, regardless of their order.
#' This is important for how terapadog handles these matching samples during
#' shuffles.
#'
#' @param df1 A dataframe.
#' @param df2 A dataframe.
#' @return NULL. Throws an error if column names do not match.
#' @keywords internal
check_matching_colnames <- function(df1, df2) {
  if (!(all(colnames(df1) %in% colnames(df2))
        && all(colnames(df2) %in% colnames(df1)))) {
    stop("SampleNames do not match between the RNA and RIBO/POLY count files.")
  }
}


#' Check the Range of Values in a Dataframe
#'
#' This function checks if the range (max - min) of numeric values in the
#' dataframe is below 1. If so, this is an indication the data was scaled
#' or anyway processed in a way that makes it not suitable for DeltaTE's
#' analysis within terapadog. The input counts must be raw counts.
#'
#' @param df A dataframe containing numeric values.
#' @return NULL. Throws an error if the range does not exceed the threshold.
#' @keywords internal
check_value_range <- function(df) {

  global_min <- min(df, na.rm = TRUE)
  global_max <- max(df, na.rm = TRUE)
  range_diff <- global_max - global_min

  # If the range is between 0 and 1, data is likely normalised/scaled.
  if (range_diff <= 1) {
    stop("Value range in the data is between 0 and 1.
         Data submitted must be RAW COUNT DATA, not normalised or scaled!")}
}

#' Check and Convert Dataframe Columns to Integers
#'
#' This function checks whether each numeric column in the provided dataframe
#' contains only integer values. If a column contains floating-point numbers,
#' it rounds the values and converts them to integers.
#'
#' @param df A dataframe to be checked and potentially modified.
#' @return A dataframe with numeric columns as integers.
#' @keywords internal
check_integer_values <- function(df) {
  are_columns_integer <- sapply(df, function(x)
  {all(x == as.integer(x), na.rm = TRUE)})
  if (!all(are_columns_integer)) {
    row_names <- rownames(df)
    df <- data.frame(lapply(df, function(x) as.integer(round(x))))
    rownames(df) <- row_names
  }
  return(df)
}

