% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_plot_dual.R
\name{plot_dual_1D}
\alias{plot_dual_1D}
\title{plot_dual_1D}
\usage{
plot_dual_1D(
  mu = (1:99)/100,
  data,
  s1,
  s2,
  type = "gauss",
  OP = FALSE,
  penalty = 2 * length(data)
)
}
\arguments{
\item{mu}{values for the dual (mu btw 0 and 1 is then re-scaled) at time n with s1 function and s2 constraint}

\item{data}{a vector of data to segment generated by \code{dataGenerator_1D}}

\item{s1}{index for q_t^s1}

\item{s2}{index for q_t^s2}

\item{type}{type of cost to use: \code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}}

\item{OP}{FALSE or TRUE, using the OP algorithm (\code{dust_R_1D}) or not to get the \code{costQ} vector, otherwise all cost set to zero}

\item{penalty}{penalty value (non-negative)}
}
\value{
plot of the dual function with threshold value in red and values leading to pruning in the dual between the green vertical lines (if there is any pruning)
}
\description{
Plotting the values of the dual (mu btw 0 and 1)
}
\examples{
data <- dataGenerator_1D(chpts = 5, parameters = 0, type = "gauss")
plot_dual_1D(data = data, s1 = 2, s2 = 1, type = "gauss")

data <- dataGenerator_1D(chpts = 10, parameters = 5, type = "poisson")
plot_dual_1D(data = data, s1 = 2, s2 = 1, type = "poisson")

data <- dataGenerator_1D(chpts = 10, parameters = 0.5, type = "geom")
plot_dual_1D(data = data, s1 = 4, s2 = 1, type = "geom")
}
