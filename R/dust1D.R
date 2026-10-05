
#' Multiple Change-Point Detection for 1D Data Using the DUST Algorithm
#'
#' Detects multiple change points in univariate time series data using the DUST algorithm.
#'
#' @param data A numeric vector representing the univariate time series. No copy of the data is made, and it is not possible to append new data for incremental analysis. For such functionality, see \code{\link{dust.object.1D}}.
#' @param penalty A finite nonnegative penalty per change point. By default,
#'   \code{2 log(length(data))}.
#' @param model A character string indicating the statistical model used for change-point detection. Default is \code{"gauss"}. Supported values include \code{"gauss"}, \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}, and \code{"variance"}.
#' @param method A character string specifying the pruning algorithm. Default is \code{"DUST"}, which is generally the most efficient in practice. The constraint index used by the pruning test is always the nearest smaller active index (the "det" rule): an earlier randomised-index variant was removed after measuring it both slower and no better at pruning. Other options include:
#' \itemize{
#'   \item \code{"DUST"}: Closed-form maximum of the decision function
#'   \item \code{"DUSTib"}: Audited one-constraint inequality test with explicit domain and boundary handling
#'   \item \code{"PELT"}: PELT pruning rule.
#'   \item \code{"OP"}: no pruning
#' }
#' @param backend Either \code{"highway"} (default) or \code{"scalar"}.
#'   \code{"highway"} uses the Highway implementation when it is available
#'   in this installation and otherwise falls back to the scalar engine.
#'   \code{"scalar"} always uses the scalar engine. The former spelling
#'   \code{"Highway"} is accepted as an alias.
#'
#' @return A list containing the information computed by the DUST algorithm.
#' \itemize{
#'   \item \code{changepoints}: the sequence of optimal change points solving our penalized optimization problem
#'   \item \code{lastIndexSet}: the last non-pruned indices at time step n (= data length)
#'   \item \code{backend}: the backend that actually ran ("highway" or "scalar") -- may differ from the requested \code{backend} if Highway was requested but is not available in this installation
#'   \item \code{nb}: vector or size n (= data length) recording the number of non-pruned indices over time
#'   \item \code{costQ}: vector or size n (= data length) recording the optimal (penalized) segmentation cost over time
#' }
#'
#' @note The input data should be first normalized by function \code{data_normalization_1D} to use the default penalty in Gaussian model, instead of value \code{2 sdDiff(data)^2 log(length(data))}.
#' The smallest index, non-pruned, is always tested for pruning with the PELT rule.
#' For every method, Binomial counts must be divided by their known number
#' of trials and Negative Binomial counts by their known size; divide the
#' penalty by the same factor to preserve the raw likelihood objective.
#' Exact-zero Exponential observations and Variance residuals are rejected
#' because their unrestricted segment likelihood is singular. Data must be
#' finite and within the chosen model's domain.
#'
#' @seealso
#' \code{\link{dataGenerator_1D}} — To generate synthetic 1D data with change points and various statistical models.
#'
#' \code{\link{data_normalization_1D}} — To normalize input data prior to applying \code{dust.1D}.
#'
#' \code{\link{dust.object.1D}} — An object-oriented version of this function that supports incremental updates via \code{append} and \code{update_partition}.
#'
#' @examples
#' ### Gaussian: 5 segments of 300 points, alternating mean shifts
#' set.seed(1)
#' true_chpts <- c(300, 600, 900, 1200, 1500)
#' y <- dataGenerator_1D(chpts = true_chpts, parameters = c(0, 2, -1, 1.5, 0.5),
#'                        sdNoise = 1, type = "gauss")
#' y <- data_normalization_1D(y, type = "gauss")
#' res <- dust.1D(data = y, model = "gauss")
#' res$changepoints  # detected change points, close to true_chpts
#'
#' ### Poisson: 4 segments of 250 points, alternating low/high rates
#' set.seed(2)
#' true_chpts <- c(250, 500, 750, 1000)
#' y <- dataGenerator_1D(chpts = true_chpts, parameters = c(2, 8, 3, 6), type = "poisson")
#' y <- data_normalization_1D(y, type = "poisson")
#' res <- dust.1D(data = y, model = "poisson")
#' res$changepoints
#'
#' ### Geometric: 3 segments of 300 points, varying success probability
#' set.seed(3)
#' true_chpts <- c(300, 600, 900)
#' y <- dataGenerator_1D(chpts = true_chpts, parameters = c(0.6, 0.2, 0.5), type = "geom")
#' y <- data_normalization_1D(y, type = "geom")
#' res <- dust.1D(data = y, model = "geom")
#' res$changepoints
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
