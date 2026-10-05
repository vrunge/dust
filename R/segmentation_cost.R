
#' Compute the Total Segmentation Cost in One Dimension with Various Models
#'
#' This function calculates the total segmentation cost for one-dimensional data using various models, such as Gaussian, Poisson, Exponential, and others.
#' It computes the segmentation cost for each segment defined by change-points (\code{chpts}) and sums these costs to provide a total segmentation cost.
#'
#' @param data A numeric vector representing the one-dimensional data to be segmented.
#' @param chpts Strictly increasing integer segment endpoints, ending at
#'   \code{length(data)}.
#' @param model A character string specifying the cost model to be used. Supported models include \code{"gauss"},
#'        \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"},
#'        and \code{"variance"}. The default model is \code{"gauss"}.
#'
#' @return A numeric value representing the total segmentation cost.
#'
#' @details
#' For each segment defined by two consecutive change-points, the function applies \code{\link{Cost_1D}}, which calculates the cost of the segment based on the provided model. Supported models include:
#'
#' \describe{
#'   \item{\code{"gauss"}}{Gaussian model, computes the negative log-likelihood under the Gaussian distribution.}
#'   \item{\code{"poisson"}}{Poisson model.}
#'   \item{\code{"exp"}}{Exponential model.}
#'   \item{\code{"geom"}}{Geometric model.}
#'   \item{\code{"bern"}}{Bernoulli model.}
#'   \item{\code{"binom"}}{Binomial model.}
#'   \item{\code{"negbin"}}{Negative Binomial model.}
#'   \item{\code{"variance"}}{Variance-based cost model.}
#' }
#'
#' For \code{"variance"}, cumulative sums of squared observations are used.
#' Binomial data must be proportions divided by the known number of trials;
#' Negative Binomial counts must be divided by the known size. For
#' \code{"gauss"}, \code{sum(data^2)/2} is added to the model's reduced cost.
#'
#' @examples
#' ### Negative Binomial series, 3 segments of 300 points with distinct
#' ### success probabilities: the total cost at the true change points is
#' ### lower than at a deliberately wrong set of change points.
#' set.seed(50)
#' true_chpts <- c(300, 600, 900)
#' data <- dataGenerator_1D(chpts = true_chpts, parameters = c(0.6, 0.2, 0.4),
#'                           nbSuccess = 10, type = "negbin")
#' segmentation_Cost_1D(data, true_chpts, model = "negbin")
#'
#' wrong_chpts <- c(150, 450, 900)
#' segmentation_Cost_1D(data, wrong_chpts, model = "negbin")  # higher cost
#'
#' @seealso \code{\link{Cost_1D}}
#'
#' @export
segmentation_Cost_1D <- function(data, chpts, model = "gauss")
{
  allowed <- c("gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance")
  if (!is.character(model) || length(model) != 1L || !model %in% allowed)
    stop("unknown model", call. = FALSE)
  if (!is.numeric(data) || !is.null(dim(data)) || length(data) == 0L || any(!is.finite(data)))
    stop("data must be a nonempty finite numeric vector", call. = FALSE)
  valid <- switch(model,
    gauss = rep(TRUE, length(data)),
    poisson = data >= 0,
    exp = data > 0,
    geom = data >= 1,
    bern = data >= 0 & data <= 1,
    binom = data >= 0 & data <= 1,
    negbin = data >= 0,
    variance = data != 0)
  if (!all(valid)) stop("data contain observations outside the model domain", call. = FALSE)
  if (!is.numeric(chpts) || length(chpts) == 0L || any(!is.finite(chpts)) ||
      any(chpts != floor(chpts)) || any(diff(chpts) <= 0L) ||
      chpts[1L] < 1L || chpts[length(chpts)] != length(data))
    stop("chpts must be strictly increasing integer endpoints ending at length(data)", call. = FALSE)
  ### we add 0 and thus move all the indices
  chpts <- c(0, chpts) + 1
  K <- length(chpts) ### K >= 2
  S <- c(0, cumsum(if (model == "variance") data^2 else data))

  totalCost <- 0
  for (i in seq.int(2L, K))
  {
    totalCost <- totalCost + Cost_1D(S, chpts[i-1], chpts[i], model)  # Apply the cost function
  }

  ### we add the sum of square in case of the Gaussian cast
  ### to get a 0 cost value in case of a no-noise perfectly well segmented data
  if(model == "gauss")
  {
    totalCost <- totalCost + sum(data^2)/2
  }
  return(totalCost)
}



#' Compute the Cost for a Single Segment Based on a Specified Model
#'
#' This function computes the cost of a single segment of data, defined by indices \code{a} and \code{b},
#' using a model specified in the \code{model} parameter.
#'
#' @param S A cumulative statistic for the data: sums of squared observations
#'   for \code{"variance"}, ordinary sums for the other models.
#' @param a An integer representing the end index of the previous segment.
#' @param b An integer representing the end index of the segment.
#' @param model A character string specifying the model to be used. Supported models include
#'        \code{"gauss"}, \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"},
#'        \code{"negbin"}, and \code{"variance"}.
#'
#' @return A numeric value representing the cost of the segment.
#'
#' @details
#' The function supports several models, including:
#'
#' \describe{
#'   \item{\code{"gauss"}}{Gaussian model, which computes the negative log-likelihood for a Gaussian distribution.}
#'   \item{\code{"poisson"}}{Poisson model.}
#'   \item{\code{"exp"}}{Exponential model.}
#'   \item{\code{"geom"}}{Geometric model.}
#'   \item{\code{"bern"}}{Bernoulli model.}
#'   \item{\code{"binom"}}{Binomial model.}
#'   \item{\code{"negbin"}}{Negative Binomial model.}
#'   \item{\code{"variance"}}{Variance-based cost model.}
#' }
#'
#' @examples
#' ### Cost of the first true segment of a 3-segment Negative Binomial series
#' set.seed(50)
#' data <- dataGenerator_1D(chpts = c(300, 600, 900), parameters = c(0.6, 0.2, 0.4),
#'                           nbSuccess = 10, type = "negbin")
#' S <- c(0, cumsum(data))
#' Cost_1D(S, a = 1, b = 301, model = "negbin")  # segment [1, 300]
#'
#' @export
Cost_1D <- function(S, a, b, model)
{
  allowed <- c("gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance")
  if (!is.character(model) || length(model) != 1L || !model %in% allowed)
    stop("unknown model", call. = FALSE)
  if (!is.numeric(S) || length(S) < 2L || any(!is.finite(S)))
    stop("S must be a finite cumulative-statistic vector", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a != floor(a) || b != floor(b) ||
      a < 1L || b > length(S) || b <= a)
    stop("a and b must define a nonempty segment within S", call. = FALSE)
  cost_value <- 0
  diff <- S[b] - S[a]
  delta <- b - a

  if(model == "gauss")
  {
    cost_value <- - 0.5 * diff^2 / delta;
  }
  if (model == "poisson")
  {
    if (diff < 0) stop("Poisson segment statistic must be nonnegative", call. = FALSE)
    if (diff != 0.0)
    {
      cost_value <- diff * (1.0 - log(diff / delta))
    }
  }
  if (model == "exp")
  {
    if (diff <= 0) stop("Exponential segment statistic must be positive", call. = FALSE)
    cost_value <- delta * (1 + log(diff / delta))
  }
  if (model == "geom")
  {
    ratio <- diff / delta
    if (ratio < 1) stop("Geometric segment mean must be at least one", call. = FALSE)
    if(ratio != 1)
    {
      cost_value <- delta * log(ratio - 1) - diff * log((ratio - 1) / ratio)
    }
  }
  if (model == "bern")
  {
    ratio <- diff / delta
    if (ratio < 0 || ratio > 1) stop("Bernoulli segment mean must be in [0,1]", call. = FALSE)
    if(ratio != 0 && ratio != 1)
    {
      cost_value <- - delta * (ratio * log(ratio) + (1 - ratio) * log(1 - ratio))
    }
  }
  if (model == "binom")
  {
    ratio <- diff / delta
    if (ratio < 0 || ratio > 1) stop("Binomial data must be proportions in [0,1]", call. = FALSE)
    if(ratio != 0 && ratio != 1)
    {
      cost_value <- - delta * (ratio * log(ratio) + (1 - ratio) * log(1 - ratio))
    }
  }
  if (model == "negbin")
  {
    ratio <- diff / delta
    if (ratio < 0) stop("Negative Binomial segment mean must be nonnegative", call. = FALSE)
    if(ratio != 0)
    {
      cost_value <- delta * log(1 + ratio) - diff * log(ratio / (1 + ratio))
    }
  }
  if (model == "variance")
  {
    if (diff <= 0) stop("Variance segment statistic must be positive", call. = FALSE)
    cost_value <- 0.5 * delta * (1.0 + log(diff / delta))
  }
  return(cost_value)
}
