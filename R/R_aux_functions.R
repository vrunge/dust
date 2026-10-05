
#########################################################
#############  backtracking_changepoint  ################
#########################################################

#' backtracking_changepoint
#'
#' @description Internal helper that walks the \code{cp} back-pointer vector
#'   produced by a segmentation algorithm to recover the ordered sequence of
#'   change points. Not meant to be called directly on raw time series data.
#' @param cp vector of changepoints of size n+1
#' @param n data length
#' @return An increasing integer vector of change-point positions, ending at \code{n}.
#' @keywords internal
backtracking_changepoint <- function(cp, n)
{
  changepoints <- n ##### vector of change-point to build
  current <- n

  while(changepoints[1] > 0)
  {
    pointval <- cp[current] ##### new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }

  changepoints <- changepoints[-1] ##### remove the first point equal to 0

  return(changepoints)
}


