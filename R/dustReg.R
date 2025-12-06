
#' Multiple Change-Point Detection for Mean and Variance Using the DUST Algorithm
#'
#' Detects multiple change points in a univariate Gaussian time series where both
#' the mean and the variance may change simultaneously, using the DUST algorithm
#' specialized to the mean/variance model.
#'
#' @param data A numeric vector representing the univariate time series.
#'   No copy of the data is made, and it is not possible to append new data for
#'   incremental analysis. For such functionality, see \code{\link{dust.object.Reg}}.
#'
#' @param penalty A positive numeric value specifying the penalty for introducing
#'   a new change point. By default, it is set to \code{4 * log(length(data))},
#'   which is a BIC-type penalty adapted to the two-parameter (mean and variance)
#'   Gaussian model.
#'
#' @param method A character string specifying the method used to handle indices
#'   and pruning tests in the mean/variance DUST algorithm. The default is
#'   \code{"det_DUST"}, which selects an efficient deterministic method for the
#'   mean/variance model.
#'
#'   In the mean/variance setting, deterministic index rules are denoted by
#'   \code{"det1"} and \code{"det2"} (replacing \code{"det"} and
#'   \code{"rand"} used in the 1D setting). The \code{method} argument controls
#'   both the index construction and the dual maximization algorithm, typically
#'   through combinations of the form \code{"det1_PRUNING"} or
#'   \code{"det2_PRUNING"}, where \code{PRUNING} indicates the pruning strategy:
#'
#' \itemize{
#'   \item \code{"DUSTr"}: Random evaluation of the dual (uniform sampling).
#'   \item \code{"DUST1"}: Inequality-based decision test using one index for
#'         pruning (internal method \code{decisionTest1}).
#'   \item \code{"DUST2"}: Improved inequality-based decision test using two
#'         indices for pruning (internal method \code{decisionTest2}).
#'   \item \code{"PELT"}: PELT pruning rule.
#'   \item \code{"OP"}: No additional DUST pruning (pure optimal partitioning).
#' }
#'
#' @return A list containing the information computed by the mean/variance DUST
#'   algorithm:
#' \itemize{
#'   \item \code{changepoints}: the sequence of optimal change points solving the
#'         penalized optimization problem;
#'   \item \code{lastIndexSet}: the last non-pruned indices at time step
#'         \code{n} (= data length);
#'   \item \code{nb}: a vector of size \code{n} recording the number of
#'         non-pruned indices over time;
#'   \item \code{costQ}: a vector of size \code{n} recording the optimal
#'         (penalized) segmentation cost over time.
#' }
#'
#' @note The smallest non-pruned index is always tested for pruning (PELT-type
#'   rule), and additional pruning is performed according to the chosen
#'   \code{method}. For incremental analysis (streaming data), use
#'   \code{\link{dust.object.Reg}} instead.
#'
#' @seealso
#' \code{\link{dataGenerator_Reg}} — To generate synthetic data with
#' simultaneous mean and variance changes.
#'
#' \code{\link{dust.object.Reg}} — Object-oriented interface for the
#'   mean/variance model that supports incremental updates via
#'   \code{append_data} and \code{update_partition}.
#'
#' \code{\link{dust.1D}} — Single-parameter version for standard 1D models.
#'
#' @examples
#' y <- c(rnorm(300, 0, 1),
#'       rnorm(400, 0, 1.5),
#'       rnorm(300, 1, 1.5))
#'
#' res <- dust.Reg(y)
#' res$changepoints
#' @export
dust.Reg <- function(
    data = data
    , penalty = 4*log(length(data))
    , method = "det_DUST"
)
{
  object <- new(DUST_Reg, method)
  object$dust(data, penalty)
  return(object$get_partition())
}
