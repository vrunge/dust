% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sdDiff}
\alias{sdDiff}
\title{Calculate Standard Deviation or MAD of Differences in a Numeric Vector}
\usage{
sdDiff(y, method)
}
\arguments{
\item{y}{A numeric vector.}

\item{method}{A character string specifying the method to use.
Options are: \code{"HALL"}, \code{"MAD"}, and \code{"SD"}.
Default is \code{"HALL"}.}
}
\value{
A numeric value representing the calculated measure of variability
  according to the specified method.
  \itemize{
    \item \code{"HALL"}: Calculates the standard deviation using a specific weighted
          difference method (HALL method).
    \item \code{"MAD"}: Returns the MAD (Median Absolute Deviation) of the differences
          between consecutive elements.
    \item \code{"SD"}: Returns the standard deviation of the differences
          between consecutive elements.
  }
}
\description{
The `sdDiff` function calculates a measure of variability (standard deviation or MAD)
of a numeric vector.
It supports three methods: "HALL", "MAD", and "SD".
}
\examples{
y <- c(rnorm(500, mean = 1), rnorm(500,mean = 2))
sdDiff(y, "HALL")
sdDiff(y, "MAD")
sdDiff(y, "SD")

}
