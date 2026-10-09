#' Multiple Change-Point Detection for Multivariate Data Using the DUST Algorithm
#'
#' @description Change-point detection in independent multivariate time series with the DUST pruning rule.
#' Each row of \code{data} is a time series, all with the same model and the same change points.
#'
#' @param data a matrix (one time series per row)
#' @param penalty the penalty for a change point. By default, \code{2 * nrow(data) * log(ncol(data))}
#' @param model the model: \code{"gauss"} (default), \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"} or \code{"variance"}
#' @param method the pruning method:
#' \itemize{
#'   \item \code{"exact"} (default): Gaussian dual maximization. Other models use numerical one-constraint maxima; with two constraints, axis maxima and unbounded directions, plus an interior critical point in dimension two. In higher dimensions the non-Gaussian test is conservative, not a general exact maximizer.
#'   \item \code{"coordinateDescent"}: maximization of the decision function, one multiplier at a time
#'   \item \code{"QN"}: maximization with a quasi-Newton algorithm (BFGS with Armijo condition)
#'   \item \code{"randomEval"}: evaluation at random points
#'   \item \code{"PELT"}: PELT pruning rule
#'   \item \code{"OP"}: no pruning
#' }
#' @param constraints number of indices used in the pruning test (the largest active indices smaller than the tested index), between 1 and \code{nrow(data)}. Default is 1. With \code{"exact"}, non-Gaussian models use at most two constraints.
#' @param nbIterations number of iterations (sweeps for \code{"coordinateDescent"}, steps for \code{"QN"}, random points for \code{"randomEval"}). By default, 1 for \code{"coordinateDescent"} and 10 otherwise.
#' @param threads number of threads for the scan of the indices. By default, all the cores for \code{"OP"}, \code{"PELT"} and \code{"PELTpar"} (many indices), 1 otherwise.
#' @param epsilon stopping rule for \code{"coordinateDescent"} and \code{"QN"} when \code{nbIterations} is \code{NULL}: the search stops when the decision function increases by less than \code{epsilon} (at most 1000 iterations)
#'
#' @return A list containing the information computed by the DUST algorithm.
#' \itemize{
#'   \item \code{changepoints}: the sequence of optimal change points
#'   \item \code{lastIndexSet}: the last non-pruned indices at time step n
#'   \item \code{nb}: number of non-pruned indices over time
#'   \item \code{costQ}: optimal penalized cost on the -2 log-likelihood scale, with segmentation-independent terms omitted
#' }
#'
#' @note The pruning is safe: an index is removed only when the decision function is positive.
#' Before each search, a bound on the decision function can stop it early (no pruning possible).
#'
#' @seealso \code{\link{dust.object.MD}}, \code{\link{dataGenerator_MD}}
#'
#' @examples
#' y <- dataGenerator_MD(chpts = c(60, 120), parameters = cbind(c(0, 2), c(0, -1)), type = "gauss")
#' dust.MD(y)$changepoints
#' dust.MD(y, method = "coordinateDescent", constraints = 2, nbIterations = 20)$changepoints
#' dust.MD(y, method = "QN", constraints = 2, epsilon = 1e-8)$changepoints
#' @export
dust.MD <- function(data,
                    penalty = 2 * nrow(data) * log(ncol(data)),
                    model = "gauss", method = "exact",
                    constraints = 1L, nbIterations = NULL, epsilon = NULL,
                    threads = .default_threads(method)) {
  if (!is.matrix(data) || !is.numeric(data) ||
      nrow(data) < 1L || ncol(data) < 1L)
    stop("data must be a nonempty numeric matrix", call. = FALSE)
  if (!is.numeric(penalty) || length(penalty) != 1L ||
      !is.finite(penalty) || penalty < 0)
    stop("penalty must be a finite nonnegative number", call. = FALSE)
  object <- dust.object.MD(model, method, constraints, nbIterations, epsilon, threads)
  object$dust(data, penalty)
}

#' dust.object.MD
#'
#' @description Constructs a DUST object for multivariate data, with methods to add data and update the segmentation
#'
#' @inheritParams dust.MD
#' @param constraints number of indices used in the pruning test (1 by default). With \code{NULL}, the number of rows of the data. With \code{"exact"}, non-Gaussian models use at most two; \code{get_info()} reports this limit.
#'
#' @details The penalty is fixed at the first call of \code{append_data} (with \code{NULL}, \code{2 * nrow * log(ncol)} of this first data matrix).
#'
#' @return A DUST object with the methods
#' \itemize{
#'   \item \code{append_data(data, penalty)}: add new data
#'   \item \code{update_partition()}: update the segmentation
#'   \item \code{get_partition()}: get the segmentation
#'   \item \code{get_info()}: get information about the object
#'   \item \code{dust(data, penalty)}: append_data, update_partition and get_partition
#' }
#'
#' @examples
#' y <- dataGenerator_MD(chpts = c(50, 100), parameters = cbind(c(0, 2), c(0, -1)), type = "gauss")
#' obj <- dust.object.MD(constraints = 2)
#' obj$append_data(y[, 1:60], 4 * log(ncol(y)))
#' obj$update_partition()
#' obj$append_data(y[, 61:100], NULL)
#' obj$update_partition()
#' obj$get_partition()$changepoints
#' @export
dust.object.MD <- function(model = "gauss", method = "exact",
                           constraints = 1L, nbIterations = NULL, epsilon = NULL,
                           threads = .default_threads(method)) {
  model <- match.arg(model, c("gauss", "poisson", "exp", "geom", "bern",
                              "binom", "negbin", "variance"))
  method <- match.arg(method, c("exact", "coordinateDescent", "QN",
                               "randomEval", "PELT", "OP"))
  if (!is.null(constraints) &&
      (!is.numeric(constraints) || length(constraints) != 1L ||
       !is.finite(constraints) || constraints != floor(constraints) ||
       constraints < 1 || constraints > .Machine$integer.max))
    stop("constraints must be a positive integer or NULL", call. = FALSE)
  if (!is.null(nbIterations) &&
      (!is.numeric(nbIterations) || length(nbIterations) != 1L ||
       !is.finite(nbIterations) || nbIterations != floor(nbIterations) ||
       nbIterations < 1 || nbIterations > .Machine$integer.max))
    stop("nbIterations must be a positive integer or NULL", call. = FALSE)
  if (!is.null(epsilon) &&
      (!is.numeric(epsilon) || length(epsilon) != 1L ||
       !is.finite(epsilon) || epsilon < 0))
    stop("epsilon must be finite, nonnegative, or NULL", call. = FALSE)
  use_epsilon <- is.null(nbIterations) && !is.null(epsilon)
  if (use_epsilon && method == "randomEval")
    stop("epsilon is not available for randomEval; set nbIterations",
         call. = FALSE)
  use_epsilon <- use_epsilon &&
    method %in% c("coordinateDescent", "QN")
  iterations <- if (!is.null(nbIterations)) nbIterations else
    if (use_epsilon) 1000L else if (method == "coordinateDescent") 1L else 10L
  new(Detector, model, method,
      if (is.null(constraints)) 0L else as.integer(constraints),
      as.integer(iterations), if (use_epsilon) as.double(epsilon) else -1, as.integer(threads))
}
