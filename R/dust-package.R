#' dust: Fast Multiple Change-Point Detection
#'
#' The \pkg{dust} package implements DUST pruning for multiple change-point
#' detection. Its \code{dust.1D} and \code{dust.MD} interfaces support eight
#' one-parameter cost models for univariate and independent multivariate data.
#' The \code{dust.meanVar} interface detects changes in both the mean and
#' variance of a Gaussian series. All three offer scalar and optional Highway
#' backends and online objects for incremental analysis.
#'
#' @name dust
#' @docType package
#' @aliases dust-package
#'
#' @useDynLib dust, .registration = TRUE
#' @import Rcpp
#' @import methods
"_PACKAGE"
