#' dust: Fast Multiple Change-Point Detection in Univariate Time Series
#'
#' The \pkg{dust} package implements DUST pruning for multiple change-point
#' detection. Its \code{dust.1D} interface supports eight one-parameter cost
#' models. The \code{dust.meanVar} interface detects changes in both the mean
#' and variance of a Gaussian series, using either one or two constraints.
#' Both interfaces offer scalar and optional Highway backends and online
#' objects for incremental analysis.
#'
#' @name dust
#' @docType package
#' @aliases dust-package
#'
#' @useDynLib dust, .registration = TRUE
#' @import Rcpp
#' @import methods
"_PACKAGE"
