% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_dust_2param.R
\name{dust_R_2param_meanVar}
\alias{dust_R_2param_meanVar}
\title{dust_R_2param_meanVar}
\usage{
dust_R_2param_meanVar(data, penalty, pruningOpt)
}
\arguments{
\item{data}{a vector of data (type \code{"meanVar"}) generated by \code{dataGenerator_meanVar}}

\item{penalty}{penalty value (non-negative)}

\item{pruningOpt}{the pruning option to use. 0 = nothing, 1 = PELT, 2 = dust, 3 = PELT + dust}
}
\description{
DUST algorithm for meanVar problem
}
