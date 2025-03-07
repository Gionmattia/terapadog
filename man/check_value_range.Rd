% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing_helpers.R
\name{check_value_range}
\alias{check_value_range}
\title{Check the Range of Values in a Dataframe}
\usage{
check_value_range(df)
}
\arguments{
\item{df}{A dataframe containing numeric values.}
}
\value{
NULL. Throws an error if the range does not exceed the threshold.
}
\description{
This function checks if the range (max - min) of numeric values in the
dataframe is below 1. If so, this is an indication the data was scaled
or anyway processed in a way that makes it not suitable for DeltaTE's
analysis within terapadog. The input counts must be raw counts.
}
\examples{
#' # Internal function, code cannot be run from here.
\dontrun{
rna_file <- system.file("extdata", "rna_counts.tsv",
package = "terapadog")
rna_df <- read.table(rna_file, header = TRUE, sep = "\t")
check_value_range(rna_df)
}
}
\keyword{internal}
