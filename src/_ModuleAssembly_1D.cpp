#include <Rcpp.h>

// --- // Models // --- //
#include "1D_1_GaussModel.h"
#include "1D_2_PoissonModel.h"
#include "1D_3_ExpModel.h"
#include "1D_4_GeomModel.h"
#include "1D_5_BernModel.h"
#include "1D_6_BinomModel.h"
#include "1D_7_NegbinModel.h"
#include "1D_8_VarianceModel.h"

using namespace Rcpp;

// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_1D *newModule1D(const std::string& model,
                     const std::string& method,
                     Nullable<int> nbLoops)
{
  ///////////////////  method separation                ///////////////////
  ///////////////////  part 1 = rand = 0, det = 1       ///////////////////
  ///////////////////  part 2 = Algo type for pruning   ///////////////////

  std::vector<std::string> method_INFO;
  size_t pos = method.find('_');  // Find the position of the underscore

  if (pos != std::string::npos)
  {
    method_INFO.push_back(method.substr(0, pos));        // First part before the underscore
    method_INFO.push_back(method.substr(pos + 1));       // Second part after the underscore
  }
  else
  {
    method_INFO.push_back(method);
    method_INFO.push_back(method);
  }

  ///////////////////  DEFAULT CHOICE  = best choice ///////////////////
  std::string dualmax_algo = "DUST"; /// exact DUST
  std::string constr_index = "det"; /// deterministic choice for constraint

  if (method_INFO[0] == "rand"){constr_index = "rand";}
  else if (method_INFO[0] == "det"){constr_index = "det";}
  else if ((method_INFO[0] == "OP") || (method_INFO[0] == "PELT")){constr_index = "-";}

  if (method_INFO[1] == "DUSTr") {dualmax_algo = "DUSTr";} //algo0
  else if (method_INFO[1] == "DUST"){dualmax_algo = "DUST";} //algo1
  else if (method_INFO[1] == "DUSTgs"){dualmax_algo = "DUSTgs";} //algo2
  else if (method_INFO[1] == "DUSTbs"){dualmax_algo = "DUSTbs";} //algo3
  else if (method_INFO[1] == "DUSTqn"){dualmax_algo = "DUSTqn";} //algo4
  else if (method_INFO[1] == "PELT"){dualmax_algo = "PELT";} //algo5
  else if (method_INFO[1] == "OP"){dualmax_algo = "OP";} //algo6

  if (model == "gauss")  return new Gauss_1D(dualmax_algo, constr_index, nbLoops);
  else if (model == "poisson") return new Poisson_1D(dualmax_algo, constr_index, nbLoops);
  else if (model == "exp") return new Exp_1D(dualmax_algo, constr_index, nbLoops);
  else if (model == "geom") return new Geom_1D(dualmax_algo, constr_index, nbLoops);
  else if (model == "bern") return new Bern_1D(dualmax_algo, constr_index, nbLoops);
  else if (model == "binom") return new Binom_1D(dualmax_algo, constr_index, nbLoops);
  else if (model == "negbin") return new Negbin_1D(dualmax_algo, constr_index, nbLoops);
  else if (model == "variance") return new Variance_1D(dualmax_algo, constr_index, nbLoops);
  else return new Gauss_1D(dualmax_algo, constr_index, nbLoops); /// DEFAULT GAUSS
  return nullptr;
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_1D to R
//'
//' @name DUST_1D
//'
//' @description
//' This module exposes the \code{DUST_1D} C++ class to R, allowing you to create
//' instances of \code{DUST_1D} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULE1D)
{
  class_<DUST_1D>("DUST_1D")

    .factory<const std::string&, const std::string&, Nullable<int>>(newModule1D)

    .method("append_data", &DUST_1D::append_data)
    .method("update_partition", &DUST_1D::update_partition)
    .method("get_partition", &DUST_1D::get_partition)
    .method("get_info", &DUST_1D::get_info)
    .method("dust", &DUST_1D::dust)
  ;
}

