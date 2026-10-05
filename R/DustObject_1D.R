
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
#' @param method A character string specifying the pruning algorithm used by
#'   the dual maximization test. The default is \code{"DUST"}, which generally
#'   selects an efficient method for the chosen model. The constraint index
#'   used by the pruning test is always the nearest smaller active index (the
#'   "det" rule): an earlier randomised-index variant was removed after
#'   measuring it both slower and no better at pruning. Currently implemented
#'   options are:
#'   \itemize{
#'     \item \code{"DUST"}: Closed-form maximizer of the decision function
#'     \item \code{"DUSTib"}: Audited one-constraint inequality test with explicit domain and boundary handling
#'     \item \code{"PELT"}: PELT pruning rule.
#'     \item \code{"OP"}: OP pruning rule.
#'   }
#' @param backend Either \code{"highway"} (default) or \code{"scalar"}.
#'   \code{"highway"} uses Highway when available and the scalar engine
#'   otherwise. \code{"scalar"} always selects the scalar engine. The former
#'   spelling \code{"Highway"} remains an accepted alias.
#' @details The first nonempty append fixes the penalty. With \code{NULL},
#'   the default is \code{2 * log(first batch size)}, which can differ from
#'   the one-shot default for the full data length. Later non-\code{NULL}
#'   penalties must equal the first penalty. Empty appends do nothing.
#'   Call \code{update_partition()} after appending and before
#'   \code{get_partition()}; requesting a partition earlier is an error.
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
#'     (parameters, internal state, and actual backend);
#'   \item \code{dust}: wrapper that runs \code{append_data},
#'     \code{update_partition} and \code{get_partition} sequentially.
#' }
#'
#' @examples
#' ### Streaming scenario: a 1200-point Gaussian series with 2 true mean
#' ### shifts (at 400 and 800) arrives in 3 uneven batches; each batch is
#' ### appended and the partition is updated without reprocessing earlier data.
#' set.seed(10)
#' true_chpts <- c(400, 800, 1200)
#' y <- dataGenerator_1D(chpts = true_chpts, parameters = c(0, 1.5, -1),
#'                        sdNoise = 1, type = "gauss")
#' y <- data_normalization_1D(y, type = "gauss")
#' penalty <- 2 * log(length(y))
#'
#' ob <- dust.object.1D(model = "gauss", method = "DUST")
#' ob$append_data(y[1:500], penalty)
#' ob$update_partition()
#' ob$get_partition()$changepoints  # based on the first 500 points only
#'
#' ob$append_data(y[501:900], penalty)
#' ob$update_partition()
#' ob$get_partition()$changepoints  # resumed, now sees 900 points
#'
#' ob$append_data(y[901:1200], penalty)
#' ob$update_partition()
#' ob$get_partition()$changepoints  # all 1200 points, close to true_chpts
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
#' @description
#' Internal Highway implementation of \code{\link{dust.object.1D}}:
#' data can be appended in batches, with
#' \code{update_partition} resuming from where the previous call left off
#' instead of reprocessing the whole series.
#'
#' Uses the Google Highway SIMD engine when available for DUST, DUSTib,
#' PELT and OP. Otherwise returns a scalar object with the same interface.
#' All methods require normalized Binomial proportions, positive Exponential
#' observations and nonzero Variance residuals. Negative Binomial counts
#' must be divided by the known size. Scale penalties accordingly.
#'
#' @param model Statistical model; see \code{\link{dust.object.1D}}.
#' @param method Pruning method; DUSTib has a native Highway kernel.
#'
#' @return
#' An object with the same methods as \code{\link{dust.object.1D}}:
#' \code{append_data}, \code{update_partition}, \code{get_partition},
#' \code{get_info} and \code{dust}.
#'
#' @examples
#' ### Same streaming scenario as dust.object.1D(), but with an Exponential
#' ### series (rate shifts at 400 and 800) and the Highway-backed object.
#' # The actual backend is available from ob$get_info()$backend.
#'
#' set.seed(11)
#' true_chpts <- c(400, 800, 1200)
#' y <- dataGenerator_1D(chpts = true_chpts, parameters = c(2, 0.5, 3), type = "exp")
#' y <- data_normalization_1D(y, type = "exp")
#' penalty <- 2 * log(length(y))
#'
#' ob <- dust.object.1D(model = "exp", method = "DUST", backend = "highway")
#' ob$append_data(y[1:450], penalty)
#' ob$update_partition()
#' ob$get_partition()$changepoints
#'
#' ob$append_data(y[451:850], penalty)
#' ob$update_partition()
#' ob$get_partition()$changepoints
#'
#' ob$append_data(y[851:1200], penalty)
#' ob$update_partition()
#' ob$get_partition()$changepoints  # close to true_chpts
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
