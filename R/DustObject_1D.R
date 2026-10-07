
## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##


Rcpp::loadModule("DUSTMODULE1D", TRUE)
Rcpp::loadModule("DUSTHWMODULE1D", TRUE)
Rcpp::loadModule("DUSTMODULEmeanVar2", TRUE)

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
#' @param method the pruning method: \code{"DUST"} (default), \code{"DUSTib"}, \code{"PELT"} or \code{"OP"} (see \code{\link{dust.1D}})
#' @param backend \code{"highway"} (default) or \code{"scalar"}
#'
#' @details The penalty is fixed at the first call of \code{append_data}. With \code{NULL}, it is \code{2 log(n)} with n the size of the first data vector.
#' Call \code{update_partition()} before \code{get_partition()}.
#'
#' @return A DUST 1D object with the methods
#' \itemize{
#'   \item \code{append_data(data, penalty)}: add new data
#'   \item \code{update_partition()}: update the segmentation with the new data
#'   \item \code{get_partition()}: get the segmentation
#'   \item \code{get_info()}: get information about the object (parameters, backend)
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
    , backend = "highway"
)
{
  backend <- .dust_backend(backend)
  if (backend == "highway")
    return(dust.object.1D.HW(model, method))
  object <- new(DUST_1D, model, method)
  return(object)
}


#' dust.object.1D.HW
#'
#' @description Highway version of \code{\link{dust.object.1D}} (scalar object if Highway is not available).
#' Use \code{dust.object.1D(..., backend = "highway")}.
#'
#' @param model the model, see \code{\link{dust.object.1D}}
#' @param method the pruning method, see \code{\link{dust.object.1D}}
#'
#' @return a DUST 1D object (same methods as \code{\link{dust.object.1D}})
#'
#' @examples
#' y <- dataGenerator_1D(chpts = c(400, 800), parameters = c(2, 0.5), type = "exp")
#' ob <- dust.object.1D(model = "exp", backend = "highway")
#' ob$dust(y, 2 * log(800))$changepoints
#' @keywords internal
dust.object.1D.HW <- function(
    model = "gauss"
    , method = "DUST"
)
{
  if (DUST.1D.HW.backend() != "highway")
  {
    return(dust.object.1D(model, method, backend = "scalar"))
  }
  object <- new(DUST_1D_HW_Obj, model, method)
  return(object)
}
