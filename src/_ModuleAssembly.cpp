#include <Rcpp.h>

// --- // Models // --- //
#include "1D_A1_GaussModel.h"
#include "1D_A2_PoissonModel.h"
#include "1D_A3_ExpModel.h"
#include "1D_A4_GeomModel.h"
#include "1D_A5_BernModel.h"
#include "1D_A6_BinomModel.h"
#include "1D_A7_NegbinModel.h"

#include "2D_DUSTmeanVar.h"
#include "2D_DUSTreg.h"

using namespace Rcpp;

// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_1D *newModule1D(const std::string& model,
                     const std::string& method,
                     Nullable<double> alpha,
                     Nullable<int> nbLoops)
{
  bool use_dual_max;
  bool random_constraint;
  if (method == "randIndex_randEval")
  {
    use_dual_max = false; /// random evaluation of the dual
    random_constraint = true;  /// random choice for the unique constraint
  }
  else if (method == "randIndex_detEval")
  {
    use_dual_max = true; /// exact evaluation of the dual
    random_constraint = true; /// random choice for the unique constraint
  }
  else
  {
    use_dual_max = true;  /// exact evaluation of the dual
    random_constraint = false;  /// choice of the closest index
  }

  if (model == "gauss")
    return new Gauss_1D(use_dual_max, random_constraint, alpha, nbLoops);
  if (model == "poisson")
    return new Poisson_1D(use_dual_max, random_constraint, alpha, nbLoops);
  if (model == "exp")
    return new Exp_1D(use_dual_max, random_constraint, alpha, nbLoops);
  if (model == "geom")
    return new Geom_1D(use_dual_max, random_constraint, alpha, nbLoops);
  if (model == "bern")
    return new Bern_1D(use_dual_max, random_constraint, alpha, nbLoops);
  if (model == "binom")
    return new Binom_1D(use_dual_max, random_constraint, alpha, nbLoops);
  if (model == "negbin")
    return new Negbin_1D(use_dual_max, random_constraint, alpha, nbLoops);
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

    .factory<const std::string&, const std::string&, Nullable<double>, Nullable<int>>(newModule1D)

    .method("init_raw", &DUST_1D::init)
    .method("compute", &DUST_1D::compute)
    .method("get_partition", &DUST_1D::get_partition)
    .method("quick_raw", &DUST_1D::quick)
  ;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


DUST_meanVar *newModuleMeanVar(const std::string& method,
                               Nullable<double> alpha,
                               Nullable<int> nbLoops)
{
  bool use_dual_max;
  bool random_constraint;
  if (method == "randIndex_randEval")
  {
    use_dual_max = false; /// random evaluation of the dual
    random_constraint = true;  /// random choice for the unique constraint
  }
  else if (method == "randIndex_detEval")
  {
    use_dual_max = true; /// exact evaluation of the dual
    random_constraint = true; /// random choice for the unique constraint
  }
  else
  {
    use_dual_max = true;  /// exact evaluation of the dual
    random_constraint = false;  /// choice of the closest index
  }

  return new DUST_meanVar(use_dual_max, random_constraint, alpha, nbLoops);
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
RCPP_MODULE(DUSTMODULEMeanVar)
{
  class_<DUST_meanVar>("DUST_meanVar")

  .factory<const std::string&, Nullable<double>, Nullable<int>>(newModuleMeanVar)

  .method("init_raw", &DUST_meanVar::init)
  .method("compute", &DUST_meanVar::compute)
  .method("get_partition", &DUST_meanVar::get_partition)
  .method("quick_raw", &DUST_meanVar::quick)
  ;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


DUST_reg *newModuleReg(const std::string& method,
                       Nullable<double> alpha,
                       Nullable<int> nbLoops)
{
  bool use_dual_max;
  bool random_constraint;
  if (method == "randIndex_randEval")
  {
    use_dual_max = false; /// random evaluation of the dual
    random_constraint = true;  /// random choice for the unique constraint
  }
  else if (method == "randIndex_detEval")
  {
    use_dual_max = true; /// exact evaluation of the dual
    random_constraint = true; /// random choice for the unique constraint
  }
  else
  {
    use_dual_max = true;  /// exact evaluation of the dual
    random_constraint = false;  /// choice of the closest index
  }

  return new DUST_reg(use_dual_max, random_constraint, alpha, nbLoops);
}


// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

//' @title MyModule: Exposing DUST_reg to R
//'
//' @name DUST_reg
//'
//' @description
//' This module exposes the \code{DUST_reg} C++ class to R, allowing you to create
//' instances of \code{DUST_reg} and call its methods directly from R.
//'
//' @export
RCPP_MODULE(DUSTMODULEreg)
{
  class_<DUST_reg>("DUST_reg")

  .factory<const std::string&, Nullable<double>, Nullable<int>>(newModuleReg)

  .method("init_raw", &DUST_reg::init)
  .method("compute", &DUST_reg::compute)
  .method("get_partition", &DUST_reg::get_partition)
  .method("quick_raw", &DUST_reg::quick)
  ;
}




