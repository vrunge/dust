
## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##


Rcpp::loadModule("DUSTMODULE1D", TRUE)

Rcpp::loadModule("FLATDUST1D", TRUE) #with reserve(n + 1);
Rcpp::loadModule("FLAT2DUST1D", TRUE) #with [n+1]

Rcpp::loadModule("FLATOP1D", TRUE)

## --------------------------------- ##
## ----///////////////////////// --- ##
## --- //    1D DUST object   // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.object.1D
#'
#' @description
#' Constructs a DUST 1D object for multiple change-point detection in univariate
#' time series.
#'
#' @param model A character string specifying the model for the data. The default
#'   is \code{"gauss"}. Available models are:
#'   \itemize{
#'     \item \code{"gauss"}: Gaussian distribution with known variance.
#'     \item \code{"poisson"}: Poisson distribution, typically for count data.
#'     \item \code{"exp"}: Exponential distribution.
#'     \item \code{"geom"}: Geometric distribution.
#'     \item \code{"bern"}: Bernoulli distribution, typically for binary data.
#'     \item \code{"binom"}: Binomial distribution, for experiments with a fixed
#'       number of trials.
#'     \item \code{"negbin"}: Negative Binomial distribution, for overdispersed
#'       count data.
#'     \item \code{"variance"}: Gaussian distribution with unknown variance and
#'       zero mean.
#'   }
#'
#' @param method A character string specifying the method used to handle indices
#'   and pruning tests in the algorithm. The default is \code{"det_DUST"},
#'   which generally selects an efficient method for the chosen model. Other
#'   available methods are:
#'   \itemize{
#'     \item \code{"rand_PRUNING"}: Random index-based methods with different dual maximization algorithms
#'     \item \code{"det_PRUNING"}: Deterministic index-based methods with different dual maximization algorithms
#'   }
#'
#'   The \code{"PRUNING"} suffix controls the algorithm used to maximize the dual
#'   function. Currently implemented options are:
#'   \itemize{
#'     \item \code{"DUSTr"}: Random evaluation of the dual (uniform sampling)
#'     \item \code{"DUSTib"}: Closed-form maximizer of the decision function. An inequality-based rule
#'     \item \code{"DUST"}: Closed-form maximizer of the dual function if exists
#'     \item \code{"DUSTgs"}: Golden-section search
#'     \item \code{"DUSTbs"}: Binary search, with an early stopping rule based on the tangent line at the current point
#'     \item \code{"DUSTqn"}: Quasi-Newton method with an Armijo line-search condition (the most efficient iterative method).
#'     \item \code{"PELT"}: PELT pruning rule.
#'     \item \code{"OP"}: OP pruning rule.
#'   }
#'
#' @param nbLoops Integer; number of iterations used in the algorithm for
#'   maximizing the dual function.
#'
#' @return
#' A DUST 1D object providing the following methods:
#' \itemize{
#'   \item \code{append_data}: add new observations to the data vector to be
#'     analysed;
#'   \item \code{update_partition}: update the optimal partition after new data
#'     have been appended;
#'   \item \code{get_partition}: retrieve the optimal partition once it has been
#'     computed;
#'   \item \code{get_info}: obtain information about the current object
#'     (parameters, internal state, etc.);
#'   \item \code{dust}: wrapper that runs \code{append_data},
#'     \code{update_partition} and \code{get_partition} sequentially.
#' }
#'
#' @examples
#' ob <- dust.object.1D()
#' ob$append_data(rnorm(100),2*log(100))
#' ob$update_partition()
#' ob$get_partition()
#' ob$append_data(rnorm(100, mean = 1),2*log(100))
#' ob$update_partition()
#' ob$get_partition()
#' ob$append_data(rnorm(100, mean = 0),2*log(100))
#' ob$update_partition()
#' ob$get_partition()
dust.object.1D <- function(
    model = "gauss"
    , method = "det_DUST"
    , nbLoops = 10
)
{
  object <- new(DUST_1D, model, method, nbLoops)
  return(object)
}

