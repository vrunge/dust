/// MD DUST: Highway state, factory and Rcpp module
/// (decision function and pruning methods in the MD_*.h files)

#include <Rcpp.h>
#include <string>
#include <vector>

#include "1D_DUSTib.h"

using namespace Rcpp;

#ifdef HAVE_HIGHWAY
namespace md_hw {

struct State
{
  int model = 0;
  std::vector<double> pos, cost;
  std::vector<std::vector<double>> sums; // one array per dimension

  void push(size_t index, double value, const double* prefix)
  {
    pos.push_back(static_cast<double>(index));
    cost.push_back(value);
    for (size_t row = 0; row < sums.size(); ++row)
      sums[row].push_back(prefix[row]);
  }

  void compact(const std::vector<unsigned int>& survivors)
  {
    if (survivors.size() == pos.size()) return;
    size_t src = 0;
    for (size_t dst = 0; dst < survivors.size(); ++dst)
    {
      while (static_cast<unsigned int>(pos[src]) != survivors[dst]) ++src;
      if (src != dst)
      {
        pos[dst] = pos[src];
        cost[dst] = cost[src];
        for (auto& row : sums) row[dst] = row[src];
      }
      ++src;
    }
    pos.resize(survivors.size());
    cost.resize(survivors.size());
    for (auto& row : sums) row.resize(survivors.size());
  }

  /// VARIANCE model only: rounding errors of the sums (see 1D_Variance.h)
  std::vector<std::vector<double>> sums_lo;

  void push_variance(size_t index, double value, const double* prefix, const double* prefix_lo)
  {
    push(index, value, prefix);
    for (size_t row = 0; row < sums_lo.size(); ++row)
      sums_lo[row].push_back(prefix_lo[row]);
  }

  void compact_variance(const std::vector<unsigned int>& survivors)
  {
    if (survivors.size() == pos.size()) return;
    size_t src = 0;
    for (size_t dst = 0; dst < survivors.size(); ++dst)
    {
      while (static_cast<unsigned int>(pos[src]) != survivors[dst]) ++src;
      if (src != dst)
      {
        pos[dst] = pos[src];
        cost[dst] = cost[src];
        for (auto& row : sums) row[dst] = row[src];
        for (auto& row : sums_lo) row[dst] = row[src];
      }
      ++src;
    }
    pos.resize(survivors.size());
    cost.resize(survivors.size());
    for (auto& row : sums) row.resize(survivors.size());
    for (auto& row : sums_lo) row.resize(survivors.size());
  }
};

} // namespace md_hw

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "MD_DUST-inl.h"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "MD_DUST-inl.h"
namespace md_hw {
HWY_EXPORT(Scan);
HWY_EXPORT(ScanVariance);
}

#endif

#include "MD_Engine.h"

namespace dust_md {

template <class Model>
MDBase* make_method(const std::string& model, const std::string& method,
                    const std::string& backend, int constraints,
                    int nb_iterations, double epsilon)
{
  if (method == "OP")
    return new MDOP<Model>(model, method, backend, constraints,
                           nb_iterations, epsilon);
  return new MDPruned<Model>(model, method, backend, constraints,
                             nb_iterations, epsilon);
}

MDBase* make_md(const std::string& model, const std::string& method,
                const std::string& backend, int constraints,
                int nb_iterations, double epsilon)
{
  if (method != "coordinateDescent" && method != "randomEval" && method != "exact" &&
      method != "QN" &&
      method != "PELT" && method != "OP") stop("invalid MD method");
  if (backend != "highway" && backend != "scalar") stop("invalid MD backend");
  if (constraints < 0) stop("constraints must be positive or NULL");
  if (nb_iterations < 1) stop("nbIterations must be a positive integer");
  if (epsilon != -1.0 && (!std::isfinite(epsilon) || epsilon < 0.0))
    stop("epsilon must be finite, nonnegative, or NULL");
  if (epsilon >= 0.0 && method == "randomEval")
    stop("epsilon is not available for randomEval");
  const std::string& effective = method;
  if (model == "gauss") return make_method<GaussPolicy>(model, effective, backend, constraints, nb_iterations, epsilon);
  if (model == "poisson") return make_method<PoissonPolicy>(model, effective, backend, constraints, nb_iterations, epsilon);
  if (model == "exp") return make_method<ExpPolicy>(model, effective, backend, constraints, nb_iterations, epsilon);
  if (model == "geom") return make_method<GeomPolicy>(model, effective, backend, constraints, nb_iterations, epsilon);
  if (model == "bern") return make_method<BernPolicy>(model, effective, backend, constraints, nb_iterations, epsilon);
  if (model == "binom") return make_method<BinomPolicy>(model, effective, backend, constraints, nb_iterations, epsilon);
  if (model == "negbin") return make_method<NegbinPolicy>(model, effective, backend, constraints, nb_iterations, epsilon);
  if (model == "variance") return make_method<VariancePolicy>(model, effective, backend, constraints, nb_iterations, epsilon);
  stop("invalid MD model");
}

} // namespace dust_md

using dust_md::MDBase;
using dust_md::make_md;

RCPP_MODULE(DUSTMODULEMD)
{
  class_<MDBase>("DUST_MD")
    .factory<const std::string&, const std::string&, const std::string&, int, int, double>(make_md)
    .method("append_data", &MDBase::append_data)
    .method("update_partition", &MDBase::update_partition)
    .method("get_partition", &MDBase::get_partition)
    .method("get_info", &MDBase::get_info)
    .method("dust", &MDBase::dust);
}
