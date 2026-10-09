
## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##


Rcpp::loadModule("DUSTMODULE", TRUE)

## --------------------------------- ##
## ----///////////////////////// --- ##
## --- //    1D DUST object   // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.object.1D
#'
#' @description Constructs a DUST 1D object for multiple change-point detection in univariate time series.
#' Data can be added step by step with \code{append_data} and \code{update_partition}.
#'
#' @param model the model for the data. Available models are:
#' \itemize{
#'   \item \code{"gauss"}: Gaussian distribution with known variance (default)
#'   \item \code{"poisson"}: Poisson distribution
#'   \item \code{"exp"}: Exponential distribution
#'   \item \code{"geom"}: Geometric distribution
#'   \item \code{"bern"}: Bernoulli distribution
#'   \item \code{"binom"}: Binomial distribution
#'   \item \code{"negbin"}: Negative Binomial distribution
#'   \item \code{"variance"}: Gaussian distribution with mean 0 and unknown variance
#' }
#' @param method the pruning method: \code{"DUST"} (default), \code{"DUSTib"}, \code{"PELT"}, \code{"PELTpar"} or \code{"OP"} (see \code{\link{dust.1D}})
#' @param threads number of threads (see \code{\link{dust.1D}})
#'
#' @details The penalty is fixed at the first call of \code{append_data}. With \code{NULL}, it is \code{2 log(n)} with n the size of the first data vector.
#' Call \code{update_partition()} before \code{get_partition()}.
#'
#' @return A DUST 1D object with the methods
#' \itemize{
#'   \item \code{append_data(data, penalty)}: add new data
#'   \item \code{update_partition()}: update the segmentation with the new data
#'   \item \code{get_partition()}: get the segmentation
#'   \item \code{get_info()}: get information about the object (parameters)
#'   \item \code{dust(data, penalty)}: append_data, update_partition and get_partition
#' }
#'
#' @examples
#' y <- dataGenerator_1D(chpts = c(400, 800, 1200), parameters = c(0, 1.5, -1), type = "gauss")
#' penalty <- 2 * log(length(y))
#'
#' ob <- dust.object.1D()
#' ob$append_data(y[1:500], penalty)
#' ob$update_partition()
#' ob$get_partition()$changepoints
#'
#' ob$append_data(y[501:1200], penalty)
#' ob$update_partition()
#' ob$get_partition()$changepoints
dust.object.1D <- function(
    model = "gauss"
    , method = "DUST"
    , threads = .default_threads(method)
)
{
  new(Detector, model, method, 1L, 1L, -1, as.integer(threads))
}
