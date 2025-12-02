
#' @export
dust.meanVar <- function(
    data = data
    , penalty = 4*log(length(data))
    , method = "det_DUST"
)
{
  object <- new(DUST_meanVar, method)
  object$dust(data, penalty)
  return(object$get_partition())
}


