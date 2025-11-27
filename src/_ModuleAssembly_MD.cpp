// --- // Models // --- //
#include "MD_A_DUST.h"
#include "MD_A1_GaussModel.h"
#include "MD_A2_PoissonModel.h"
#include "MD_A3_ExpModel.h"
#include "MD_A4_GeomModel.h"
#include "MD_A5_BernModel.h"
#include "MD_A6_BinomModel.h"
#include "MD_A7_NegbinModel.h"
#include "MD_A8_VarianceModel.h"

using namespace Rcpp;

// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_MD *newModuleMD(const std::string& model,
                     const std::string& method,
                     Nullable<unsigned> nbLoops)
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
  std::string dualmax_algo = "DUSTqn"; /// quasi newton
  std::string constr_index = "det"; /// deterministic choice

  if (method_INFO[0] == "rand"){constr_index = "rand";} // rand index choice
  else if (method_INFO[0] == "det"){constr_index = "det";} // det index choice

  if (method_INFO[1] == "DUSTr") {dualmax_algo = "DUSTr";} //algo0
  else if (method_INFO[1] == "DUST"){dualmax_algo = "DUST";} //algo1
  else if (method_INFO[1] == "DUSTgs"){dualmax_algo = "DUSTgs";} //algo2
  else if (method_INFO[1] == "DUSTbs"){dualmax_algo = "DUSTbs";} //algo3
  else if (method_INFO[1] == "DUSTqn"){dualmax_algo = "DUSTqn";} //algo4
  else if (method_INFO[1] == "PELT"){dualmax_algo = "PELT";} //algo5
  else if (method_INFO[1] == "OP"){dualmax_algo = "OP";} //algo6

  if (model == "gauss") return new Gauss_MD(dualmax_algo, constr_index, nbLoops);
  else if (model == "poisson") return new Poisson_MD(dualmax_algo, constr_index, nbLoops);
  else if (model == "exp") return new Exp_MD(dualmax_algo, constr_index, nbLoops);
  else if (model == "geom") return new Geom_MD(dualmax_algo, constr_index, nbLoops);
  else if (model == "bern") return new Bern_MD(dualmax_algo, constr_index, nbLoops);
  else if (model == "binom") return new Binom_MD(dualmax_algo, constr_index, nbLoops);
  else if (model == "negbin") return new Negbin_MD(dualmax_algo, constr_index, nbLoops);
  else if (model == "variance") return new Variance_MD(dualmax_algo, constr_index, nbLoops);
  else return new Gauss_MD(dualmax_algo, constr_index, nbLoops); /// DEFAULT GAUSS
  return nullptr;
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_MD to R
//'
//' @name DUST_MD
//'
//' @description
//' This module exposes the \code{DUST_MD} C++ class to R, allowing you to create
//' instances of \code{DUST_MD} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULEMD)
{
  class_<DUST_MD>("DUST_MD")

  .factory<const std::string&, const std::string&, Nullable<unsigned>>(newModuleMD)

  .method("append_data", &DUST_MD::append_data)
  .method("update_partition", &DUST_MD::update_partition)
  .method("get_partition", &DUST_MD::get_partition)
  .method("get_info", &DUST_MD::get_info)
  .method("dust", &DUST_MD::dust)
  ;
}




