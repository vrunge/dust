Rcpp::loadModule("DUSTMODULEMD", TRUE)

#' Detect Changes in Independent Multivariate Data
#'
#' Each column of \code{data} is one observation and each row is an independent
#' component. All components use the same cost model. Segment costs are added
#' across rows, while change points are shared.
#'
#' @param data A nonempty numeric matrix with components in rows and time in
#'   columns.
#' @param penalty A finite nonnegative penalty per change point. By default,
#'   \code{2 * nrow(data) * log(ncol(data))}.
#' @param model One of \code{"gauss"}, \code{"poisson"}, \code{"exp"},
#'   \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}, or
#'   \code{"variance"}; see \code{\link{dust.1D}}.
#' @param method \code{"coordinateDescent"} (default), \code{"iterative"},
#'   \code{"QN"}, \code{"randomEval"}, \code{"exact"}, \code{"PELT"}, or
#'   \code{"OP"}. Coordinate descent minimizes the negative of the concave
#'   DUST decision function. \code{"iterative"} uses projected gradient ascent
#'   with backtracking. \code{"QN"} uses safeguarded inverse BFGS updates and
#'   an Armijo line search. Both search jointly over the selected constraints
#'   and support all eight models, checking the mean domain at every trial.
#'   Multipliers forced to zero by a boundary segment mean are held fixed.
#'   Random evaluation checks sampled feasible points. For \code{model = "gauss"},
#'   \code{"exact"} solves the joint decision problem over all selected
#'   constraints by checking stationary faces of a concave quadratic.
#'   It tries a short active-face search, then enumerates faces if needed.
#'   The fallback can require exponentially many faces as \code{constraints}
#'   grows. If numerical rank or optimality checks are inconclusive, the candidate
#'   is retained. For other models, \code{"exact"} uses \code{"PELT"}.
#'   These methods prune only
#'   when a feasible decision value is strictly positive with a numerical
#'   tolerance. A failed search or exhausted budget retains the candidate;
#'   a finite search is not guaranteed to find the maximum. \code{"PELT"}
#'   uses only the decision at zero, and \code{"OP"} does no pruning.
#' @param backend \code{"highway"} (default) or \code{"scalar"}. Highway is used if
#'   available in the installed package; otherwise the scalar engine runs.
#' @param constraints Number of earlier active indices used for each DUST
#'   test, from 1 to \code{nrow(data)}. The default uses all \code{nrow(data)}.
#'   They are always the largest active indices smaller than the tested
#'   index. All search methods use the available earlier indices, up to
#'   this limit, even when fewer than \code{constraints} exist.
#' @param nbIterations Positive integer or \code{NULL} (default). A separate
#'   budget applies to each candidate's pruning test at each time point.
#'   For \code{"coordinateDescent"}, one iteration is one sweep over all
#'   selected multiplier coordinates. For \code{"iterative"}, it is one
#'   projected-gradient step with up to 60 backtracking trials. For
#'   \code{"QN"}, it is one quasi-Newton step with up to 60 backtracking
#'   trials and, if needed, one gradient fallback with up to 60 more.
#'   For \code{"randomEval"}, it is one random draw. Ignored by
#'   \code{"exact"}, \code{"PELT"}, and \code{"OP"}. When \code{NULL}, the
#'   effective budget is 10 unless \code{epsilon} is active.
#'   An explicit \code{nbIterations} takes priority over \code{epsilon}.
#' @param epsilon \code{NULL} (default) or a finite nonnegative threshold
#'   for the absolute gain in the normalized decision function after a
#'   complete coordinate sweep or an accepted optimizer step. Used only when
#'   \code{nbIterations}
#'   is \code{NULL}, with a cap of 1000 sweeps or iterations per candidate.
#'   The search stops when the gain is at most \code{epsilon}; this is a
#'   stopping heuristic, not a certificate that the maximum was found.
#'   Supported by \code{"coordinateDescent"}, \code{"iterative"}, and
#'   \code{"QN"}. Other methods ignore it, except \code{"randomEval"},
#'   which rejects an active \code{epsilon} because random draws have no
#'   meaningful consecutive gain.
#'
#' @return A list with \code{changepoints}, \code{lastIndexSet}, actual \code{backend},
#'   \code{nb} (active candidates over time), and \code{costQ} (optimal costs).
#' @seealso \code{\link{dust.object.MD}}
#' @examples
#' set.seed(13)
#' y <- rbind(c(rnorm(60), rnorm(60, 2)),
#'            c(rnorm(60), rnorm(60, -1)))
#' dust.MD(y, model = "gauss", constraints = 2)$changepoints
#' dust.MD(y, method = "iterative", constraints = 2, nbIterations = 20)$changepoints
#' dust.MD(y, method = "QN", constraints = 2, epsilon = 1e-8)$changepoints
#' @export
dust.MD <- function(data,
                    penalty = 2 * nrow(data) * log(ncol(data)),
                    model = "gauss", method = "coordinateDescent",
                    backend = "highway", constraints = nrow(data),
                    nbIterations = NULL, epsilon = NULL) {
  if (!is.matrix(data) || !is.numeric(data) ||
      nrow(data) < 1L || ncol(data) < 1L)
    stop("data must be a nonempty numeric matrix", call. = FALSE)
  if (!is.numeric(penalty) || length(penalty) != 1L ||
      !is.finite(penalty) || penalty < 0)
    stop("penalty must be a finite nonnegative number", call. = FALSE)
  object <- dust.object.MD(model, method, backend, constraints,
                           nbIterations, epsilon)
  object$dust(data, penalty)
}

#' Create an Incremental Multivariate DUST Object
#'
#' Append matrix batches with a fixed number of rows. Call
#' \code{update_partition()} after appending and before \code{get_partition()}.
#' The first nonempty append fixes the penalty. If it is \code{NULL}, the default
#' is \code{2 * nrow(first batch) * log(ncol(first batch))}; later non-\code{NULL}
#' penalties must equal it.
#'
#' @inheritParams dust.MD
#' @param constraints Number of earlier active indices. \code{NULL} (default)
#'   uses the number of rows in the first nonempty batch. With Gaussian
#'   \code{method = "exact"}, these indices are tested jointly.
#' @return An object with \code{append_data(data, penalty)}, \code{update_partition()},
#'   \code{get_partition()}, \code{get_info()}, and \code{dust(data, penalty)} methods.
#' @examples
#' set.seed(14)
#' y <- rbind(c(rnorm(50), rnorm(50, 2)),
#'            c(rnorm(50), rnorm(50, -1)))
#' obj <- dust.object.MD(constraints = 2)
#' obj$append_data(y[, 1:60], 4 * log(ncol(y)))
#' obj$update_partition()
#' obj$append_data(y[, 61:100], NULL)
#' obj$update_partition()
#' obj$get_partition()$changepoints
#' @export
dust.object.MD <- function(model = "gauss", method = "coordinateDescent",
                           backend = "highway", constraints = NULL,
                           nbIterations = NULL, epsilon = NULL) {
  model <- match.arg(model, c("gauss", "poisson", "exp", "geom", "bern",
                              "binom", "negbin", "variance"))
  method <- match.arg(method, c("coordinateDescent", "iterative", "QN",
                               "randomEval", "exact", "PELT", "OP"))
  backend <- .dust_backend(backend)
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
    method %in% c("coordinateDescent", "iterative", "QN")
  iterations <- if (!is.null(nbIterations)) nbIterations else
    if (use_epsilon) 1000L else 10L
  new(DUST_MD, model, method, backend,
      if (is.null(constraints)) 0L else as.integer(constraints),
      as.integer(iterations), if (use_epsilon) as.double(epsilon) else -1)
}
