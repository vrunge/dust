
#' Multiple Change-Point Detection for 1D Data Using the DUST Algorithm
#'
#' Detects multiple change points in univariate time series data using the DUST algorithm.
#'
#' @param data A numeric vector representing the univariate time series. No copy of the data is made, and it is not possible to append new data for incremental analysis. For such functionality, see \code{\link{dust.object.1D}}.
#' @param penalty A positive numeric value specifying the penalty for introducing a new change point. By default, it is set to \code{2 log(length(data))}.
#' @param model A character string indicating the statistical model used for change-point detection. Default is \code{"gauss"}. Supported values include \code{"gauss"}, \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}, and \code{"variance"}.
#' @param method A character string specifying the index-handling and pruning strategies. Default is \code{"detIndex_Eval4"}, which selects an efficient deterministic method for the chosen model. Other options include:
#'
#' \itemize{
#'   \item \code{"rand_PRUNING"}: Random index selection methods that randomly choose an index smaller than the tested index \code{s}, using different dual maximization algorithms
#'   \item \code{"det_PRUNING"}: Deterministic index selection methods that pick the largest index smaller than tested index \code{s}, also using different dual maximization algorithms
#' }
#'
#' Dual maximization algorithms (\code{PRUNING}) used in the method names are defined as follows. \code{DUST} is generally the most efficient in practice. At each iteration, the smallest non-pruned index is tested using the PELT pruning rule:
#' \itemize{
#'   \item \code{"DUSTr"}: Random evaluation using a uniform distribution.
#'   \item \code{"DUSTib"}: Closed-form maximizer of the decision function. An inequality-based rule
#'   \item \code{"DUST"}: Closed-form maximum (for Gaussian model only). Otherwise, no pruning is applied, leading to the slower OP algorithm.
#'   \item \code{"DUSTgs"}: Golden-section search.
#'   \item \code{"DUSTbs"}: Binary search with early stopping based on tangent evaluation.
#'   \item \code{"DUSTqn"}: Quasi-Newton method with Armijo condition.
#'   \item \code{"PELT"}: PELT pruning rule.
#'   \item \code{"OP"}: no pruning
#' }
#' @param nbLoops An integer. The number of loops to run in the max dual optimization algorithm. Default is 10.
#'
#'' @return A list containing the information computed by the DUST algorithm.
#' \itemize{
#'   \item \code{changepoints}: the sequence of optimal change points solving our penalized optimization problem
#'   \item \code{lastIndexSet}: the last non-pruned indices at time step n (= data length)
#'   \item \code{nb}: vector or size n (= data length) recording the number of non-pruned indices over time
#'   \item \code{costQ}: vector or size n (= data length) recording the optimal (penalized) segmentation cost over time
#' }
#'
#' @note The input data should be first normalized by function \code{data_normalization_1D} to use the default penalty in Gaussian model, instead of value \code{2 sdDiff(data)^2 log(length(data))}.
#' The smallest index, non-pruned, is always tested for pruning with the PELT rule.
#'
#' @seealso
#' \code{\link{dataGenerator_1D}} — To generate synthetic 1D data with change points and various statistical models.
#'
#' \code{\link{data_normalization_1D}} — To normalize input data prior to applying \code{dust.1D}.
#'
#' \code{\link{dust.object.1D}} — An object-oriented version of this function that supports incremental updates via \code{append} and \code{update_partition}.
#'
#' @examples
#' y <- rnorm(100)
#' y <- data_normalization_1D(y)
#' dust.1D(y)
#'
#' y <- dataGenerator_1D(chpts = c(50,100), parameters = c(0.5,0.2), type = "geom")
#' y <- data_normalization_1D(y, type = "geom")
#' dust.1D(data = y, model = "geom")
#'
#' y <- dataGenerator_1D(chpts = c(20,70), parameters = c(20,5), type = "poisson")
#' y <- data_normalization_1D(y, type = "poisson")
#' dust.1D(data = y, model = "poisson")
#' @export
dust.1D <- function(
    data = data
    , penalty = 2*log(length(data))
    , model = "gauss"
    , method = "det_DUST"
    , nbLoops = 10
)
{
  object <- new(DUST_1D, model, method, nbLoops)
  object$dust(data, penalty)
  return(object$get_partition())
}


