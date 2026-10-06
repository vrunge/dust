#include <Rcpp.h>

#include "1D_DUST_Impl.h"
#include "1D_OP_Impl.h"
#include "1D_Indices.h"
#include "1D_DualMaxPolicies.h"
#include "1D_1_GaussModel.h"
#include "1D_2_PoissonModel.h"
#include "1D_3_ExpModel.h"
#include "1D_4_GeomModel.h"
#include "1D_5_BernModel.h"
#include "1D_6_BinomModel.h"
#include "1D_7_NegbinModel.h"
#include "1D_8_VarianceModel.h"

using namespace Rcpp;

template <class Model, class DualMaxPolicy>
DUST_1D *newModuleT()
{
  return new DUST_1D_T<Model, DualMaxPolicy>();
}

template <class Model>
DUST_1D *newModuleByAlgo(const std::string& method)
{
  if (method == "DUST")   return newModuleT<Model, DualMax_DUST<Model>>();
  if (method == "DUSTib") return newModuleT<Model, DualMax_DUSTib<Model>>();
  if (method == "PELT")   return newModuleT<Model, DualMax_PELT<Model>>();
  // OP scans all previous endpoints in its own engine.
  if (method == "OP")     return new DUST_1D_OP_T<Model>();
  stop("Unrecognized method \"" + method + "\". Valid values: \"DUST\", \"DUSTib\", \"PELT\", \"OP\".");
}

DUST_1D *newModule1D(const std::string& model,
                     const std::string& method)
{
  if (model == "gauss")  return newModuleByAlgo<GaussPolicy>(method);
  else if (model == "poisson") return newModuleByAlgo<PoissonPolicy>(method);
  else if (model == "exp") return newModuleByAlgo<ExpPolicy>(method);
  else if (model == "geom") return newModuleByAlgo<GeomPolicy>(method);
  else if (model == "bern") return newModuleByAlgo<BernPolicy>(method);
  else if (model == "binom") return newModuleByAlgo<BinomPolicy>(method);
  else if (model == "negbin") return newModuleByAlgo<NegbinPolicy>(method);
  else if (model == "variance") return newModuleByAlgo<VariancePolicy>(method);
  else stop("Unrecognized model \"" + model + "\".");
}


//' @title Scalar segmentation module
//'
//' @name DUST_1D
//'
//' @description
//' Rcpp module for the scalar DUST segmentation object.
//'
//' @export
RCPP_MODULE(DUSTMODULE1D)
{
  class_<DUST_1D>("DUST_1D")

    .factory<const std::string&, const std::string&>(newModule1D)

    .method("append_data", &DUST_1D::append_data)
    .method("update_partition", &DUST_1D::update_partition)
    .method("get_partition", &DUST_1D::get_partition)
    .method("get_info", &DUST_1D::get_info)
    .method("dust", &DUST_1D::dust)
  ;
}
