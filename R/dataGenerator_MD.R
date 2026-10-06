
####################################################
#############    dataGenerator_MD   ################
####################################################

#' dataGenerator_MD
#'
#' @description Generating independent multivariate time series with shared
#'   change points using the models of \code{\link{dataGenerator_1D}}.
#' @param chpts A vector of strictly increasing positive integer change-point
#'   indices. The last value is the data length.
#' @param parameters A numeric matrix or data frame with one row per segment
#'   and one column per component. For \code{"variance"}, parameters are
#'   standard deviations; for \code{"exp"}, they are rates.
#' @param sdNoise (type \code{"gauss"}) Noise standard deviation, either one
#'   value shared by all components or one value per component.
#' @param gamma (type \code{"gauss"}) Exponential decay factors in \eqn{(0,1]}.
#'   A scalar applies to all segments and components. A vector of length
#'   \code{length(chpts)} gives segment factors shared by all components.
#'   A numeric matrix or data frame has one column per component and either
#'   one row or one row per segment. Decay restarts at each segment.
#'   The default, 1, gives piecewise constant means.
#' @param nbTrials (type \code{"binom"}) Positive integer number of trials,
#'   either one value shared by all components or one value per component.
#' @param nbSuccess (type \code{"negbin"}) Positive integer number of successes,
#'   either one value shared by all components or one value per component.
#' @param type The common model: \code{"gauss"}, \code{"exp"},
#'   \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"},
#'   \code{"negbin"}, or \code{"variance"}.
#' @return A numeric matrix with components in rows and observations in
#'   columns, the orientation expected by \code{\link{dust.MD}}. Its dimensions are
#'   \code{ncol(parameters)} by the last value of \code{chpts}.
#'   A single component still returns a matrix. Components are generated
#'   independently by successive calls to \code{\link{dataGenerator_1D}}.
#' @details As in \code{dataGenerator_1D}, the binomial and negative-binomial
#'   models return raw counts. Before fitting \code{dust.MD}, normalize these
#'   componentwise with \code{\link{data_normalization_1D}}, using the
#'   corresponding \code{nbTrials} or \code{nbSuccess} as \code{size}.
#' @examples
#' dataGenerator_MD(chpts = c(50, 100),
#'                  parameters = data.frame(ts1 = c(0, 1), ts2 = c(2, -1)),
#'                  sdNoise = c(0.2, 0.5), type = "gauss")
#' dataGenerator_MD(chpts = c(50, 100),
#'                  parameters = data.frame(ts1 = c(10, 10), ts2 = c(5, 20)),
#'                  gamma = data.frame(ts1 = c(1, 0.9), ts2 = c(0.8, 1)),
#'                  sdNoise = 0.2, type = "gauss")
#' dataGenerator_MD(chpts = c(50, 100),
#'                  parameters = cbind(c(3, 5), c(2, 7)), type = "poisson")
#' dataGenerator_MD(chpts = c(50, 100),
#'                  parameters = cbind(c(0.4, 0.7), c(0.8, 0.2)),
#'                  nbTrials = c(10, 30), type = "binom")
#' dataGenerator_MD(chpts = c(50, 100),
#'                  parameters = cbind(c(1, 3), c(2, 0.5)), type = "variance")
#' @export
dataGenerator_MD <- function(chpts = 100,
                             parameters = data.frame(ts1 = 0.5, ts2 = 0.5),
                             sdNoise = 1,
                             gamma = 1,
                             nbTrials = 10,
                             nbSuccess = 10,
                             type = "gauss")
{
  if(!is.numeric(chpts) || length(chpts) == 0L ||
     any(!is.finite(chpts)) || any(chpts <= 0) || any(chpts != floor(chpts)))
    stop('chpts must be finite positive integers')
  if(is.unsorted(chpts, strictly = TRUE))
    stop('chpts should be a strictly increasing vector of change-point positions (indices)')

  numeric_matrix <- function(value, name)
  {
    if(!(is.matrix(value) || is.data.frame(value)))
      stop(name, ' must be a numeric matrix or data frame')
    if(is.data.frame(value) && !all(vapply(value, is.numeric, logical(1))))
      stop(name, ' values must all be numeric')
    value <- as.matrix(value)
    if(!is.numeric(value) || any(!is.finite(value)))
      stop(name, ' values must all be finite and numeric')
    value
  }

  parameters <- numeric_matrix(parameters, 'parameters')
  p <- ncol(parameters)
  if(p == 0L || nrow(parameters) != length(chpts))
    stop('parameters must have one row per segment and at least one column')

  allowed.types <- c("gauss", "exp", "poisson", "geom", "bern", "binom", "negbin", "variance")
  if(!is.character(type) || length(type) != 1L || is.na(type) || !type %in% allowed.types)
    stop('type must be one of: ', paste(allowed.types, collapse = ", "))

  per_component <- function(value, name)
  {
    if(!is.numeric(value) || !is.null(dim(value)) || any(!is.finite(value)) ||
       !length(value) %in% c(1L, p))
      stop(name, ' must be a finite numeric vector of length 1 or the number of components')
    rep(value, length.out = p)
  }

  if(type == "gauss")
  {
    sdNoise <- per_component(sdNoise, 'sdNoise')
    if(is.matrix(gamma) || is.data.frame(gamma))
    {
      gamma <- numeric_matrix(gamma, 'gamma')
      if(ncol(gamma) != p || !nrow(gamma) %in% c(1L, length(chpts)))
        stop('gamma must have one column per component and one row or one row per segment')
    }
    else
    {
      if(!is.numeric(gamma) || !is.null(dim(gamma)) || any(!is.finite(gamma)) ||
         !length(gamma) %in% c(1L, length(chpts)))
        stop('gamma must be a finite numeric vector of length 1 or the length of chpts')
      gamma <- matrix(gamma, nrow = length(gamma), ncol = p)
    }
  }
  if(type == "binom"){nbTrials <- per_component(nbTrials, 'nbTrials')}
  if(type == "negbin"){nbSuccess <- per_component(nbSuccess, 'nbSuccess')}

  res <- matrix(NA_real_, nrow = p, ncol = chpts[length(chpts)])
  for(i in seq_len(p))
  {
    args <- list(chpts = chpts, parameters = parameters[, i], type = type)
    if(type == "gauss")
    {
      args$sdNoise <- sdNoise[i]
      args$gamma <- gamma[, i]
    }
    if(type == "binom"){args$nbTrials <- nbTrials[i]}
    if(type == "negbin"){args$nbSuccess <- nbSuccess[i]}
    res[i, ] <- do.call(dataGenerator_1D, args)
  }
  return(res)
}
