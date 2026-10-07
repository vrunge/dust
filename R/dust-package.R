#' dust: Fast Multiple Change-Point Detection
#'
#' Multiple change-point detection with the DUST pruning rule (DUality Simple Test).
#' \code{dust.1D} for univariate data (8 models), \code{dust.MD} for independent multivariate data
#' and \code{dust.meanVar} for changes in mean and variance. Each function has an object version
#' to add data step by step (\code{dust.object.1D}, \code{dust.object.MD}, \code{dust.object.meanVar}).
#'
#' @name dust
#' @docType package
#' @aliases dust-package
#'
#' @useDynLib dust, .registration = TRUE
#' @import Rcpp
#' @import methods
"_PACKAGE"
