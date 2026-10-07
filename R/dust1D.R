
#' Multiple Change-Point Detection for 1D Data Using the DUST Algorithm
#'
#' @description Detects multiple change points in univariate time series data using the DUST algorithm.
#'
#' @param data a numeric vector (the univariate time series). To add data later, see \code{\link{dust.object.1D}}.
#' @param penalty the penalty for a change point. By default, \code{2 log(length(data))}.
#' @param model the model: \code{"gauss"} (default), \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"} or \code{"variance"}
#' @param method the pruning method. Default is \code{"DUST"}, which is generally the most efficient in practice.
#' \itemize{
#'   \item \code{"DUST"}: closed-form maximum of the decision function
#'   \item \code{"DUSTib"}: one-constraint test with explicit domain checks
#'   \item \code{"PELT"}: PELT pruning rule
#'   \item \code{"OP"}: no pruning
#' }
#' The constraint is always the largest non-pruned index smaller than the tested index.
#' @param backend \code{"highway"} (default) or \code{"scalar"}. The scalar engine is used if Highway is not available.
#'
#' @return A list containing the information computed by the DUST algorithm.
#' \itemize{
#'   \item \code{changepoints}: the sequence of optimal change points solving our penalized optimization problem
#'   \item \code{lastIndexSet}: the last non-pruned indices at time step n (= data length)
#'   \item \code{backend}: the backend used ("highway" or "scalar")
#'   \item \code{nb}: vector of size n (= data length) recording the number of non-pruned indices over time
#'   \item \code{costQ}: vector of size n (= data length) recording the optimal (penalized) segmentation cost over time
#' }
#'
#' @note The input data should be first normalized by function \code{data_normalization_1D} to use the default penalty in Gaussian model, instead of value \code{2 sdDiff(data)^2 log(length(data))}.
#' The smallest index, non-pruned, is always tested for pruning with the PELT rule.
#' Binomial and negative binomial data have to be divided by the number of trials (or successes): see \code{data_normalization_1D}.
#'
#' @seealso
#' \code{\link{dataGenerator_1D}} to generate data,
#' \code{\link{data_normalization_1D}} to normalize the data,
#' \code{\link{dust.object.1D}} to add data step by step.
#'
#' @examples
#' y <- dataGenerator_1D(chpts = c(300, 600, 900), parameters = c(0, 2, -1), type = "gauss")
#' y <- data_normalization_1D(y)
#' dust.1D(y)$changepoints
#'
#' y <- dataGenerator_1D(chpts = c(250, 500, 750), parameters = c(2, 8, 3), type = "poisson")
#' y <- data_normalization_1D(y, type = "poisson")
#' dust.1D(y, model = "poisson")$changepoints
#'
#' y <- dataGenerator_1D(chpts = c(300, 600), parameters = c(0.6, 0.2), type = "geom")
#' dust.1D(y, model = "geom")$changepoints
#' @export
dust.1D <- function(
    data
    , penalty = 2*log(length(data))
    , model = "gauss"
    , method = "DUST"
    , backend = "highway"
)
{
  backend <- .dust_backend(backend)
  if (!is.numeric(data) || !is.null(dim(data)) || length(data) == 0L)
    stop("data must be a nonempty numeric vector", call. = FALSE)
  if (!is.numeric(penalty) || length(penalty) != 1L || !is.finite(penalty) || penalty < 0)
    stop("penalty must be a finite nonnegative number", call. = FALSE)
  if (backend == "highway")
    return(DUST.1D.HW(data, penalty, model, method))
  object <- new(DUST_1D, model, method)
  object$dust(data, penalty)
  return(object$get_partition())
}
