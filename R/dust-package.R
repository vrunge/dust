#' dust: Multiple Change-Point Detection in Multivariate Time Series
#'
#' The \pkg{dust} package implements algorithms for detecting multiple change
#' points in multivariate time series by minimizing a penalized likelihood.
#' It relies on optimal partitioning dynamic programming combined with the
#' DUality Simple Test (DUST) pruning method.
#'
#' @name dust
#' @docType package
#' @aliases dust-package
#'
#' @useDynLib dust, .registration = TRUE
#' @import Rcpp
#' @import methods
"_PACKAGE"
