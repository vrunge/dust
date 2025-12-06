
## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULEmeanVar", TRUE)


## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // mean/var DUST object // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.object.meanVar
#'
#' @description
#' Constructs a DUST object for multiple change-point detection in a Gaussian
#' model with (simultaneous) changes in mean and variance ("meanVar" model).
#'
#' @param method A character string specifying the method used to handle indices
#'   and pruning tests in the algorithm. The default is \code{"det2_DUST2"},
#'   which generally selects an efficient method for the mean/variance model.
#'   Other available methods follow the pattern \code{"detK_PRUNING"}, where
#'   \code{K} indicates the index construction rule and
#'   \code{PRUNING} controls the dual maximization algorithm:
#'   \itemize{
#'     \item \code{"det1_PRUNING"}: Deterministic index-based method with
#'       constraints of type 1 (internally \code{Indices_2D_Det1});
#'     \item \code{"det2_PRUNING"}: Deterministic index-based method with
#'       constraints of type 2 (internally \code{Indices_2D_Det2}).
#'   }
#'
#'   The \code{"PRUNING"} part controls the algorithm used to maximize the dual
#'   function. Currently implemented options in the mean/variance setting are:
#'   \itemize{
#'     \item \code{"DUSTr"}: Random evaluation of the dual (uniform sampling);
#'     \item \code{"DUST1"}: Inequality-based decision test using one index
#'       for pruning (method \code{decisionTest1});
#'     \item \code{"DUST2"}: Improved inequality-based decision test using two
#'       indices for pruning (method \code{decisionTest2});
#'     \item \code{"PELT"}: PELT pruning rule;
#'     \item \code{"OP"}: OP pruning rule (no additional DUST pruning).
#'   }
#'
#'   Typical choices are for example \code{"det2_DUST2"} (default),
#'   \code{"det1_DUST1"}, \code{"det2_PELT"}, \code{"det1_DUSTr"}, etc.
#'
#' @return
#' A \code{DUST_meanVar} object providing at least the methods:
#' \itemize{
#'   \item \code{append_data(y, penalty)}: add new observations to the data
#'     vector (numeric). If \code{penalty} is \code{NULL}, the C++ code uses
#'     \code{2 * log(n)} at first call, with \code{n} the total length;
#'   \item \code{update_partition()}: run the dynamic programming and pruning
#'     updates after new data have been appended;
#'   \item \code{get_partition()}: returns a list with fields
#'     \code{changepoints}, \code{lastIndexSet}, \code{nb}, \code{costQ};
#'   \item \code{get_info()}: internal statistics (cumsums, penalty, etc.);
#'   \item \code{dust(y, penalty)}: one-shot wrapper that runs
#'     \code{append_data}, \code{update_partition} and \code{get_partition}.
#' }
#'
#' @examples
#' ob <- dust.object.meanVar()
#' ob$append_data( c(rnorm(50, 0, 1), rnorm(50, 1, 3), rnorm(50, 0, 1)), 2 * log(150))
#' ob$update_partition()
#' part <- ob$get_partition()
#' part$changepoints
#' ob$append_data( c(rnorm(50, 0, 1), rnorm(50, 1, 3), rnorm(50, 0, 1)), 2 * log(150))
#' ob$update_partition()
#' part <- ob$get_partition()
#' part$changepoints
dust.object.meanVar <- function(method = "det2_DUST2")
{
  object <- new(DUST_meanVar, method)
  return(object)
}



################################################################################

## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULEReg", TRUE)


## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // mean/var DUST object // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.object.meanVar
#'
#' @description
#' Constructs a DUST object for multiple change-point detection in a Gaussian
#' model with (simultaneous) changes in mean and variance ("meanVar" model).
#'
#' @param method A character string specifying the method used to handle indices
#'   and pruning tests in the algorithm. The default is \code{"det2_DUST2"},
#'   which generally selects an efficient method for the mean/variance model.
#'   Other available methods follow the pattern \code{"detK_PRUNING"}, where
#'   \code{K} indicates the index construction rule and
#'   \code{PRUNING} controls the dual maximization algorithm:
#'   \itemize{
#'     \item \code{"det1_PRUNING"}: Deterministic index-based method with
#'       constraints of type 1 (internally \code{Indices_2D_Det1});
#'     \item \code{"det2_PRUNING"}: Deterministic index-based method with
#'       constraints of type 2 (internally \code{Indices_2D_Det2}).
#'   }
#'
#'   The \code{"PRUNING"} part controls the algorithm used to maximize the dual
#'   function. Currently implemented options in the mean/variance setting are:
#'   \itemize{
#'     \item \code{"DUSTr"}: Random evaluation of the dual (uniform sampling);
#'     \item \code{"DUST1"}: Inequality-based decision test using one index
#'       for pruning (method \code{decisionTest1});
#'     \item \code{"DUST2"}: Improved inequality-based decision test using two
#'       indices for pruning (method \code{decisionTest2});
#'     \item \code{"PELT"}: PELT pruning rule;
#'     \item \code{"OP"}: OP pruning rule (no additional DUST pruning).
#'   }
#'
#'   Typical choices are for example \code{"det2_DUST2"} (default),
#'   \code{"det1_DUST1"}, \code{"det2_PELT"}, \code{"det1_DUSTr"}, etc.
#'
#' @return
#' A \code{DUST_meanVar} object providing at least the methods:
#' \itemize{
#'   \item \code{append_data(y, penalty)}: add new observations to the data
#'     vector (numeric). If \code{penalty} is \code{NULL}, the C++ code uses
#'     \code{2 * log(n)} at first call, with \code{n} the total length;
#'   \item \code{update_partition()}: run the dynamic programming and pruning
#'     updates after new data have been appended;
#'   \item \code{get_partition()}: returns a list with fields
#'     \code{changepoints}, \code{lastIndexSet}, \code{nb}, \code{costQ};
#'   \item \code{get_info()}: internal statistics (cumsums, penalty, etc.);
#'   \item \code{dust(y, penalty)}: one-shot wrapper that runs
#'     \code{append_data}, \code{update_partition} and \code{get_partition}.
#' }
#'
#' @examples
## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULEReg", TRUE)


## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // mean/var DUST object // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.object.Reg
#'
#' @description
#' Constructs a DUST object for multiple change-point detection in a Gaussian
#' model with (simultaneous) changes in mean and variance ("Reg" model).
#'
#' @param method A character string specifying the method used to handle indices
#'   and pruning tests in the algorithm. The default is \code{"det2_DUST2"},
#'   which generally selects an efficient method for the mean/variance model.
#'   Other available methods follow the pattern \code{"detK_PRUNING"}, where
#'   \code{K} indicates the index construction rule and
#'   \code{PRUNING} controls the dual maximization algorithm:
#'   \itemize{
#'     \item \code{"det1_PRUNING"}: Deterministic index-based method with
#'       constraints of type 1 (internally \code{Indices_2D_Det1});
#'     \item \code{"det2_PRUNING"}: Deterministic index-based method with
#'       constraints of type 2 (internally \code{Indices_2D_Det2}).
#'   }
#'
#'   The \code{"PRUNING"} part controls the algorithm used to maximize the dual
#'   function. Currently implemented options in the mean/variance setting are:
#'   \itemize{
#'     \item \code{"DUSTr"}: Random evaluation of the dual (uniform sampling);
#'     \item \code{"DUST1"}: Inequality-based decision test using one index
#'       for pruning (method \code{decisionTest1});
#'     \item \code{"DUST2"}: Improved inequality-based decision test using two
#'       indices for pruning (method \code{decisionTest2});
#'     \item \code{"PELT"}: PELT pruning rule;
#'     \item \code{"OP"}: OP pruning rule (no additional DUST pruning).
#'   }
#'
#'   Typical choices are for example \code{"det2_DUST2"} (default),
#'   \code{"det1_DUST1"}, \code{"det2_PELT"}, \code{"det1_DUSTr"}, etc.
#'
#' @return
#' A \code{DUST_Reg} object providing at least the methods:
#' \itemize{
#'   \item \code{append_data(y, penalty)}: add new observations to the data
#'     vector (numeric). If \code{penalty} is \code{NULL}, the C++ code uses
#'     \code{2 * log(n)} at first call, with \code{n} the total length;
#'   \item \code{update_partition()}: run the dynamic programming and pruning
#'     updates after new data have been appended;
#'   \item \code{get_partition()}: returns a list with fields
#'     \code{changepoints}, \code{lastIndexSet}, \code{nb}, \code{costQ};
#'   \item \code{get_info()}: internal statistics (cumsums, penalty, etc.);
#'   \item \code{dust(y, penalty)}: one-shot wrapper that runs
#'     \code{append_data}, \code{update_partition} and \code{get_partition}.
#' }
#'
#' @examples
#' ob <- dust.object.Reg()
#' ob$append_data( c(rnorm(50, 0, 1), rnorm(50, 1, 3), rnorm(50, 0, 1)), 2 * log(150))
#' ob$update_partition()
#' part <- ob$get_partition()
#' part$changepoints
#' ob$append_data( c(rnorm(50, 0, 1), rnorm(50, 1, 3), rnorm(50, 0, 1)), 2 * log(150))
#' ob$update_partition()
#' part <- ob$get_partition()
#' part$changepoints
dust.object.Reg <- function(method = "det2_DUST2")
{
  object <- new(DUST_Reg, method)
  return(object)
}
#' ob <- dust.object.Reg()
#' ob$append_data( c(rnorm(50, 0, 1), rnorm(50, 1, 3), rnorm(50, 0, 1)), 2 * log(150))
#' ob$update_partition()
#' part <- ob$get_partition()
#' part$changepoints
#' ob$append_data( c(rnorm(50, 0, 1), rnorm(50, 1, 3), rnorm(50, 0, 1)), 2 * log(150))
#' ob$update_partition()
#' part <- ob$get_partition()
#' part$changepoints
dust.object.Reg <- function(method = "det2_DUST2")
{
  object <- new(DUST_Reg, method)
  return(object)
}
