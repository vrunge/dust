#' Detect Changes in Both Gaussian Mean and Variance
#'
#' @param data A finite numeric vector.
#' @param penalty Nonnegative penalty per change point. The default is
#'   \code{4 * log(length(data))}.
#' @param method \code{"1D"} (default) tests one constraint; \code{"2D"} tests two.
#'   Both methods fit the same mean-and-variance Gaussian model.
#' @param backend \code{"highway"} (default) uses SIMD when available and otherwise
#'   falls back to \code{"scalar"}; \code{"scalar"} always uses the scalar engine.
#'   The former spelling \code{"Highway"} remains an accepted alias.
#' @return A list with \code{changepoints}, \code{lastIndexSet}, actual \code{backend},
#'   \code{nb} (active candidates by time), and \code{costQ} (optimal costs by time).
#' @seealso \code{\link{dust.object.meanVar}}
#' @examples
#' set.seed(1)
#' y <- c(rnorm(80), rnorm(80, mean = 2, sd = 1.5))
#' dust.meanVar(y, method = "1D")$changepoints
#' dust.meanVar(y, method = "2D")$changepoints
#' @export
dust.meanVar <- function(data, penalty = 4 * log(length(data)),
                         method = "1D", backend = "highway") {
  method <- match.arg(method, c("1D", "2D"))
  backend <- .dust_backend(backend)
  if (!is.numeric(data) || !is.null(dim(data)) || length(data) == 0L || any(!is.finite(data)))
    stop("data must be a nonempty finite numeric vector", call. = FALSE)
  if (!is.numeric(penalty) || length(penalty) != 1L || !is.finite(penalty) || penalty < 0)
    stop("penalty must be a finite nonnegative number", call. = FALSE)
  object <- dust.object.meanVar(method, backend)
  object$dust(data, penalty)
}

#' Create an Online Mean-and-Variance DUST Object
#'
#' The object has \code{append_data(data, penalty)}, \code{update_partition()},
#' \code{get_partition()}, \code{get_info()}, and \code{dust(data, penalty)} methods.
#' Append batches and call \code{update_partition()} to resume the computation.
#' The penalty is fixed at the first nonempty append. If its \code{penalty} argument
#' is \code{NULL}, that penalty is \code{4 * log(first batch size)}; later supplied
#' penalties must equal the first penalty. The model needs at least two
#' nonidentical values in a segment to have finite cost; an object with no
#' finite complete segmentation reports an error from \code{get_partition()}.
#' Call \code{update_partition()} after each append before retrieving a partition.
#'
#' @inheritParams dust.meanVar
#' @return An online DUST object.
#' @examples
#' set.seed(2)
#' y <- c(rnorm(80), rnorm(80, mean = 2, sd = 1.5))
#' ob <- dust.object.meanVar(method = "2D")
#' ob$append_data(y[1:80], 4 * log(length(y)))
#' ob$update_partition()
#' ob$append_data(y[81:160], NULL)
#' ob$update_partition()
#' ob$get_partition()$changepoints
#' @export
dust.object.meanVar <- function(method = "1D", backend = "highway") {
  method <- match.arg(method, c("1D", "2D"))
  backend <- .dust_backend(backend)
  new(DUST_meanVar2, method, if (backend == "highway") "highway" else "scalar")
}
