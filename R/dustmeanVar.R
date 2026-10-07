#' Change-Point Detection in Mean and Variance Using the DUST Algorithm
#'
#' @description Detection of changes in mean and variance in a Gaussian time series with the DUST pruning rule (two-parameter model).
#'
#' @param data a numeric vector
#' @param penalty the penalty for a change point. By default, \code{4 log(length(data))}
#' @param method \code{"1D"} (default, one constraint), \code{"2D"} (two constraints) or \code{"PELT"}
#' @param backend \code{"highway"} (default) or \code{"scalar"}
#'
#' @return A list with \code{changepoints}, \code{lastIndexSet}, \code{backend}, \code{nb} (number of non-pruned indices over time) and \code{costQ} (optimal costs over time)
#'
#' @note A segment needs at least two different values to have a finite cost.
#'
#' @seealso \code{\link{dust.object.meanVar}}
#'
#' @examples
#' y <- c(rnorm(80), rnorm(80, mean = 2, sd = 1.5))
#' dust.meanVar(y)$changepoints
#' dust.meanVar(y, method = "2D")$changepoints
#' @export
dust.meanVar <- function(data, penalty = 4 * log(length(data)),
                         method = "1D", backend = "highway") {
  method <- match.arg(method, c("1D", "2D", "PELT"))
  backend <- .dust_backend(backend)
  if (!is.numeric(data) || !is.null(dim(data)) || length(data) == 0L || any(!is.finite(data)))
    stop("data must be a nonempty finite numeric vector", call. = FALSE)
  if (!is.numeric(penalty) || length(penalty) != 1L || !is.finite(penalty) || penalty < 0)
    stop("penalty must be a finite nonnegative number", call. = FALSE)
  object <- dust.object.meanVar(method, backend)
  object$dust(data, penalty)
}

#' dust.object.meanVar
#'
#' @description Constructs a DUST object for changes in mean and variance, with methods to add data and update the segmentation
#'
#' @inheritParams dust.meanVar
#'
#' @details The penalty is fixed at the first call of \code{append_data} (with \code{NULL}, \code{4 log(n)} with n the size of the first data).
#'
#' @return A DUST object with the methods \code{append_data(data, penalty)}, \code{update_partition()}, \code{get_partition()}, \code{get_info()} and \code{dust(data, penalty)}
#'
#' @examples
#' y <- c(rnorm(80), rnorm(80, mean = 2, sd = 1.5))
#' ob <- dust.object.meanVar(method = "2D")
#' ob$append_data(y[1:80], 4 * log(length(y)))
#' ob$update_partition()
#' ob$append_data(y[81:160], NULL)
#' ob$update_partition()
#' ob$get_partition()$changepoints
#' @export
dust.object.meanVar <- function(method = "1D", backend = "highway") {
  method <- match.arg(method, c("1D", "2D", "PELT"))
  backend <- .dust_backend(backend)
  new(DUST_meanVar2, method, if (backend == "highway") "highway" else "scalar")
}
