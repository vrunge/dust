
Rcpp::loadModule("DUSTMODULEMeanVar", TRUE)

## ---------------------------- ##
## ----//////////////////// --- ##
## --- // 2D DUST Object // --- ##
## ----//////////////////// --- ##
## ---------------------------- ##

#' dust.object.meanVar
#'
#' @description Generates a DUST object
#' @param method the method for handling the indices and pruning tests in the algorithm. defaults to "detIndex_Eval4", which returns the fastest method available on input model. available methods: "randIndex_randEval" prunes indices based on a random dual test, where the value of the dual is evaluated at a random point, and the dual is defined by a random index; "randIndex_detEval", where the dual is defined by a random index but evaluated at its maximum.
#' @param nbLoops number of iteration in the optimization algorithm for maximizing the dual function
#' @return a DUST object object that provides methods : fit, for fitting the data; compute, once fit has been called, for computing the optimal partition of the data; get_partition, for retrieving the optimal partition once it has been computed; and quick, a wrapper for the other 3 methods.
#' @examples
#' dust.object.meanVar()
dust.object.meanVar <- function(
    method = "detIndex_Eval4"
    , nbLoops = 10
)
{
  object <- new(DUST_meanVar, method, nbLoops)

  assign(
    "init",
    function(data, penalty = NULL)
      object$init_raw(data, penalty),
    envir = object
  )

  assign(
    "quick",
    function(data, penalty = NULL)
      object$quick_raw(data, penalty),
    envir = object
  )

  return(object)
}


####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULEreg", TRUE)


## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // 2D DUST object // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.object.reg
#'
#' @description Generates a DUST object
#' @param method the method for handling the indices and pruning tests in the algorithm.
#' defaults to "detIndex_Eval4", which returns the fastest method available on input model.
#' available methods: "randIndex_randEval" prunes indices based on a random dual test,
#' where the value of the dual is evaluated at a random point, and the dual is defined
#' by a random index; "randIndex_detEval", where the dual is defined by a random index
#' but evaluated at its maximum.
#' @param nbLoops number of iteration in the optimization algorithm for maximizing the dual function
#' @return a DUST object object that provides methods : fit, for fitting the data; compute, once fit has been called, for computing the optimal partition of the data; get_partition, for retrieving the optimal partition once it has been computed; and quick, a wrapper for the other 3 methods.
#' @examples
#' dust.object.reg()
dust.object.reg <- function(
    method = "detIndex_Eval4"
    , nbLoops = 10
)
{
  object <- new(DUST_reg, method, nbLoops)

  assign(
    "init",
    function(data, penalty = NULL)
      object$init_raw(data, penalty),
    envir = object
  )

  assign(
    "quick",
    function(data, penalty = NULL)
      object$quick_raw(data, penalty),
    envir = object
  )

  return(object)
}
