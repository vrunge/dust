
#' Multiple Change-Point Detection for 1D Data Using the DUST Algorithm
#'
#' @description Detects multiple change points in univariate time series data using the DUST algorithm.
#'
#' @param data a numeric vector (the univariate time series). To add data later, see \code{\link{dust.object.1D}}.
#' @param penalty the penalty for a change point. By default, \code{2 log(length(data))} (BIC: one parameter per segment and the change location).
#' @param model the model: \code{"gauss"} (default), \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"} or \code{"variance"}
#' @param method the pruning method. Default is \code{"DUST"}, which is generally the most efficient in practice.
#' \itemize{
#'   \item \code{"DUST"}: closed-form maximum of the decision function
#'   \item \code{"DUSTib"}: one-constraint test with explicit domain checks
#'   \item \code{"PELT"}: PELT pruning rule
#'   \item \code{"OP"}: no pruning
#'   \item \code{"PELTpar"}: PELT rule tested only on the smallest candidates, so they stay contiguous (faster scan, with threads)
#' }
#' The constraint is always the largest non-pruned index smaller than the tested index.
#' @param threads number of threads for the scan of the indices. By default, all the cores for \code{"OP"}, \code{"PELT"} and \code{"PELTpar"} (many indices), 1 otherwise.
#' @param size number of trials (\code{"binom"}) or of successes (\code{"negbin"}), \code{NULL} otherwise. With \code{size}, \code{data} are the raw counts: they are divided by \code{size} and the penalty is divided by \code{size} too, so that the default penalty and \code{costQ} are on the -2 log-likelihood scale. With \code{NULL}, the data have to be the proportions (binomial) or the counts divided by the number of successes (negative binomial), and the penalty is used as is.
#'
#' @return A list containing the information computed by the DUST algorithm.
#' \itemize{
#'   \item \code{changepoints}: the sequence of optimal change points solving our penalized optimization problem
#'   \item \code{lastIndexSet}: the last non-pruned indices at time step n (= data length)
#'   \item \code{nb}: vector of size n (= data length) recording the number of non-pruned indices over time
#'   \item \code{costQ}: vector of size n recording the optimal penalized cost on the -2 log-likelihood scale, with segmentation-independent terms omitted. For Gaussian data, add \code{sum(data[1:t]^2)} to \code{costQ[t]} to obtain the residual sum of squares plus penalties.
#' }
#'
#' @note The input data should be first normalized by function \code{data_normalization_1D} to use the default penalty in Gaussian model, instead of value \code{2 sdDiff(data)^2 log(length(data))}.
#' The smallest index, non-pruned, is always tested for pruning with the PELT rule.
#' For binomial and negative binomial data, give the raw counts and \code{size}; without \code{size}, divide the data by the number of trials (or successes), see \code{data_normalization_1D}, and the penalty too.
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
#' dust.1D(y, model = "poisson")$changepoints
#'
#' y <- dataGenerator_1D(chpts = c(300, 600), parameters = c(0.6, 0.2), type = "geom")
#' dust.1D(y, model = "geom")$changepoints
#'
#' y <- dataGenerator_1D(chpts = c(300, 600), parameters = c(0.6, 0.2), nbTrials = 10, type = "binom")
#' dust.1D(y, model = "binom", size = 10)$changepoints
#' @export
dust.1D <- function(
    data
    , penalty = 2*log(length(data))
    , model = "gauss"
    , method = "DUST"
    , threads = .default_threads(method)
    , size = NULL
)
{
  if (!is.numeric(data) || !is.null(dim(data)) || length(data) == 0L)
    stop("data must be a nonempty numeric vector", call. = FALSE)
  if (!is.numeric(penalty) || length(penalty) != 1L || !is.finite(penalty) || penalty < 0)
    stop("penalty must be a finite nonnegative number", call. = FALSE)
  dust.object.1D(model, method, threads, size)$dust(data, penalty)
}
