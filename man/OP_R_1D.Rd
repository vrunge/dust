% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_OP.R
\name{OP_R_1D}
\alias{OP_R_1D}
\title{OP_R_1D}
\usage{
OP_R_1D(data, penalty = 2 * log(length(data)), type = "gauss")
}
\arguments{
\item{data}{a vector of data to segment generated by \code{dataGenerator_1D}}

\item{penalty}{penalty value (non-negative)}

\item{type}{type of cost to use: \code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}}
}
\value{
a list with the change-point elements (each last index of each segment) and a vector \code{costQ} saving the optimal cost at each time step
}
\description{
OP algorithm for univariate time-series (with different possible data models)
}
\examples{
 OP_R_1D(dataGenerator_1D(chpts = c(200,400), c(0,1), sdNoise = 1, type = "gauss"), log(400))
 OP_R_1D(dataGenerator_1D(chpts = c(50,100), c(1,7), type = "exp"), log(100))
 OP_R_1D(dataGenerator_1D(chpts = c(50,100), c(3,10), type = "poisson"), 10*log(100))
 OP_R_1D(dataGenerator_1D(chpts = c(50,100), c(0.7,0.3), type = "geom"), 5*log(100))
 OP_R_1D(dataGenerator_1D(chpts = c(50,100), c(0.7,0.2), type = "bern"), log(100)/2)
 OP_R_1D(dataGenerator_1D(chpts = c(50,100), c(0.7, 0.3), nbTrials = 5, type = "binom"), 5*log(100))
 OP_R_1D(dataGenerator_1D(chpts = c(50,100), c(0.4,0.7),
         nbSuccess = 10, type = "negbin"), 50*log(100))
}
