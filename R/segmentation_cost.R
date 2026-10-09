
#' segmentation_Cost_1D
#'
#' @description Total cost of a segmentation (sum of the segment costs computed by \code{\link{Cost_1D}})
#'
#' @param data a numeric vector
#' @param chpts the change points (increasing, the last one is \code{length(data)})
#' @param model the model: \code{"gauss"} (default), \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}, \code{"variance"}
#'
#' @return the cost of the segmentation (without penalty)
#'
#' @details Costs are -2 log-likelihoods (the scale of the default penalty). For \code{"gauss"}, \code{sum(data^2)} is added (up to a constant).
#' Binomial and negative binomial data have to be divided by the number of trials (or successes).
#'
#' @examples
#' data <- dataGenerator_1D(chpts = c(300, 600, 900), parameters = c(0.6, 0.2, 0.4),
#'                          nbSuccess = 10, type = "negbin")
#' segmentation_Cost_1D(data, c(300, 600, 900), model = "negbin")
#' segmentation_Cost_1D(data, c(150, 450, 900), model = "negbin")  # higher cost
#'
#' @seealso \code{\link{Cost_1D}}
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
    if (model == "variance")
    {
      # sum of y^2 in the segment (more accurate than the cumsum on all data)
      local <- c(0, cumsum(data[chpts[i-1]:(chpts[i] - 1)]^2))
      totalCost <- totalCost + Cost_1D(local, 1, length(local), model)
    }
    else totalCost <- totalCost + Cost_1D(S, chpts[i-1], chpts[i], model)  # Apply the cost function
  }

  ### we add the sum of square in case of the Gaussian cast
  ### to get a 0 cost value in case of a no-noise perfectly well segmented data
  if(model == "gauss")
  {
    totalCost <- totalCost + sum(data^2)
  }
  return(totalCost)
}



#' Cost_1D
#'
#' @description Cost of the segment (a, b] from the cumulative sums of the data
#'
#' @param S the cumulative sums (of the data, or of the squared data for \code{"variance"}) with \code{S[1] = 0}
#' @param a index of the beginning of the segment in \code{S}
#' @param b index of the end of the segment in \code{S}
#' @param model the model: \code{"gauss"}, \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}, \code{"variance"}
#'
#' @return the cost of the segment
#'
#' @examples
#' data <- dataGenerator_1D(chpts = c(300, 600), parameters = c(0.6, 0.2),
#'                          nbSuccess = 10, type = "negbin")
#' S <- c(0, cumsum(data))
#' Cost_1D(S, a = 1, b = 301, model = "negbin")  # data[1:300]
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
  return(2 * cost_value)   # -2 log-likelihood
}
