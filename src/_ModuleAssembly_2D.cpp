#include <Rcpp.h>

// --- // Models // --- //
#include "2DmeanVar_DUST.h"
#include "2DReg_DUST.h"

using namespace Rcpp;


// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_meanVar *newModulemeanVar(const std::string& method)
{
  ///////////////////  method separation                ///////////////////
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
  std::string dualmax_algo = "DUST1"; /// exact DUST
  std::string constr_index = "det1"; /// deterministic choice for constraint

  if (method_INFO[0] == "det1"){constr_index = "det1";}
  else if (method_INFO[0] == "det2"){constr_index = "det2";}
  else if ((method_INFO[0] == "OP") || (method_INFO[0] == "PELT")){constr_index = "-";}

  if (method_INFO[1] == "DUSTr") {dualmax_algo = "DUSTr";} //algo0
  else if (method_INFO[1] == "DUST1"){dualmax_algo = "DUST1";} //decision test 1D
  else if (method_INFO[1] == "DUST2"){dualmax_algo = "DUST2";} //decision test 1D
  else if (method_INFO[1] == "PELT"){dualmax_algo = "PELT";} //algo1
  else if (method_INFO[1] == "OP"){dualmax_algo = "OP";} //algo1

  return new DUST_meanVar(dualmax_algo, constr_index);
  return nullptr;
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_meanVar to R
//'
//' @name DUST_meanVar
//'
//' @description
//' This module exposes the \code{DUST_meanVar} C++ class to R, allowing you to create
//' instances of \code{DUST_meanVar} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULEmeanVar)
{
  class_<DUST_meanVar>("DUST_meanVar")

    .factory<const std::string&>(newModulemeanVar)

    .method("append_data", &DUST_meanVar::append_data)
    .method("update_partition", &DUST_meanVar::update_partition)
    .method("get_partition", &DUST_meanVar::get_partition)
    .method("get_info", &DUST_meanVar::get_info)
    .method("dust", &DUST_meanVar::dust)
  ;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




DUST_Reg *newModuleReg(const std::string& method)
{
  ///////////////////  method separation                ///////////////////
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
  std::string dualmax_algo = "DUST1"; /// exact DUST
  std::string constr_index = "det1"; /// deterministic choice for constraint

  if (method_INFO[0] == "det1"){constr_index = "det1";}
  else if (method_INFO[0] == "det2"){constr_index = "det2";}
  else if ((method_INFO[0] == "OP") || (method_INFO[0] == "PELT")){constr_index = "-";}

  if (method_INFO[1] == "DUSTr") {dualmax_algo = "DUSTr";} //algo0
  else if (method_INFO[1] == "DUST1"){dualmax_algo = "DUST1";} //decision test 1D
  else if (method_INFO[1] == "DUST2"){dualmax_algo = "DUST2";} //decision test 1D
  else if (method_INFO[1] == "PELT"){dualmax_algo = "PELT";} //algo1
  else if (method_INFO[1] == "OP"){dualmax_algo = "OP";} //algo1

  return new DUST_Reg(dualmax_algo, constr_index);
  return nullptr;
}

// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_Reg to R
//'
//' @name DUST_Reg
//'
//' @description
//' This module exposes the \code{DUST_Reg} C++ class to R, allowing you to create
//' instances of \code{DUST_Reg} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULEReg)
{
   class_<DUST_Reg>("DUST_Reg")

   .factory<const std::string&>(newModuleReg)

   .method("append_data", &DUST_Reg::append_data)
   .method("update_partition", &DUST_Reg::update_partition)
   .method("get_partition", &DUST_Reg::get_partition)
   .method("get_info", &DUST_Reg::get_info)
   .method("dust", &DUST_Reg::dust)
   ;
}
