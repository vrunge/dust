% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dust1D.R
\name{dust.1D}
\alias{dust.1D}
\title{Run the 1D Dust Change Point Detection}
\usage{
dust.1D(
  data = data,
  penalty = 2 * log(length(data)),
  model = "gauss",
  method = "fastest",
  alpha = 1e-09,
  nbLoops = 10
)
}
\arguments{
\item{data}{A numeric vector. The time series data on which change point detection is performed.}

\item{penalty}{A numeric value. The penalty applied for adding a new change point. By default, it is set to \code{2 * log(length(data))}.}

\item{model}{A character string. Specifies the model used for change point detection. Default is "gauss". Possible values could include "gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance".}

\item{method}{A character string specifying the method used to handle indices and pruning tests in the algorithm. The default is \code{fastest}, which automatically selects the quickest method for the chosen model. Other available methods are:
\itemize{
  \item \code{"randIndex_Eval0"} to \code{"randIndex_Eval6"}: Random index-based methods with different dual maximization algorithm (0 through 5).
  \item \code{"detIndex_Eval0"} to \code{"detIndex_Eval6"}: Deterministic index-based methods  with different dual maximization algorithm (0 through 5).
}
Here are the current available algorithms (\code{Eval4} is often the most efficient one)
\itemize{
  \item \code{"Eval0"}: random evaluation of the dual (with uniform distribution)
  \item \code{"Eval1"}: max value with closed formula (gauss model only), otherwise no pruning performed and we get the (slow) OP algorithm
  \item \code{"Eval2"}: golden-section search.
  \item \code{"Eval3"}: binary search. At each step, we evaluate the tangent line to the current point at its max to stop the search at early step (when possible)
  \item \code{"Eval4"}: quasi-Newton method with armijo condition
  \item \code{"Eval5"}: PELT rule
  \item \code{"Eval6"}: OP rule
}}

\item{alpha}{A numeric value. For randomness level. Default is \code{1e-9}.}

\item{nbLoops}{An integer. The number of loops to run in the max dual optimization algorithm. Default is 10.}
}
\value{
A list containing the change points detected by the DUST algorithm.
}
\description{
This function performs change point detection on one-dimensional data using the DUST algorithm.
}
\examples{
data <- rnorm(1000)
result <- dust.1D(data = data)
}
