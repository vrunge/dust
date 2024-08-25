#include <Rcpp.h>

// --- // Models // --- //
#include "1D_A1_GaussModel.h"
#include "1D_A2_PoissonModel.h"
#include "1D_A3_ExpModel.h"
#include "1D_A4_GeomModel.h"
#include "1D_A5_BernModel.h"
#include "1D_A6_BinomModel.h"
#include "1D_A7_NegbinModel.h"

using namespace Rcpp;


// ---------------------------- //
// --- //////////////////// --- //
// --- // Object factory // --- //
// --- //////////////////// --- //
// ---------------------------- //

DUST_1D *newModule1D(const std::string& model, const std::string& method, Nullable<double> alpha)
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
    return new Gauss_1D(use_dual_max, random_constraint, alpha);
  if (model == "poisson")
    return new Poisson_1D(use_dual_max, random_constraint, alpha);
  if (model == "exp")
    return new Exp_1D(use_dual_max, random_constraint, alpha);
  if (model == "geom")
    return new Geom_1D(use_dual_max, random_constraint, alpha);
  if (model == "bern")
    return new Bern_1D(use_dual_max, random_constraint, alpha);
  if (model == "binom")
    return new Binom_1D(use_dual_max, random_constraint, alpha);
  if (model == "negbin")
    return new Negbin_1D(use_dual_max, random_constraint, alpha);
  return nullptr;
}



// --------------------------------- //
// --- ///////////////////////// --- //
// --- // Exposing the module // --- //
// --- ///////////////////////// --- //
// --------------------------------- //

RCPP_MODULE(DUSTMODULE1D)
{
  class_<DUST_1D>("DUST_1D")

    .factory<const std::string&, const std::string&, Nullable<double>>(newModule1D)

    .method("init_raw", &DUST_1D::init)
    .method("compute", &DUST_1D::compute)
    .method("get_partition", &DUST_1D::get_partition)
    .method("quick_raw", &DUST_1D::quick)
  ;
}