/// 1D DUST with Highway (SIMD)
/// candidates in contiguous arrays, scalar engine if no Highway

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <memory>
#include <tuple>
#include <type_traits>

using namespace Rcpp;

#include "1D_DUST.h"

#ifdef HAVE_HIGHWAY
#include "1D_DualMaxPolicies.h"
#include "1D_Validation.h"

namespace dust_hw {

enum Method { DUST, PELT, OP, DUSTIB };

// same order as dustib::model_id
using Models = std::tuple<GaussPolicy, PoissonPolicy, ExpPolicy, GeomPolicy,
                          BernPolicy, BinomPolicy, NegbinPolicy, VariancePolicy>;
template <int K> using ScalarModel = std::tuple_element_t<K, Models>;

// f(std::integral_constant<int, K>{}) for the model id K
template <class F>
void with_model(int K, F&& f)
{
  switch (K)
  {
    case 0: return f(std::integral_constant<int, 0>{});
    case 1: return f(std::integral_constant<int, 1>{});
    case 2: return f(std::integral_constant<int, 2>{});
    case 3: return f(std::integral_constant<int, 3>{});
    case 4: return f(std::integral_constant<int, 4>{});
    case 5: return f(std::integral_constant<int, 5>{});
    case 6: return f(std::integral_constant<int, 6>{});
    case 7: return f(std::integral_constant<int, 7>{});
  }
}

////////////////////////////////////////////////////////////////////////////////
/// state of the dynamic programming
/// candidate s = A_pos[i]: sum A_cumsum[i], cost A_cost[i],
/// mean statistic A_b[i] and mean cost A_q[i] since its predecessor
struct State
{
  int model = -1, method = -1;
  double penalty = 0.0;
  unsigned n = 0;
  std::vector<double> cumsum, costRecord;
  std::vector<int> chptRecord, nb;
  std::vector<double> A_cumsum, A_cost, A_pos, A_b, A_q, val;
  std::vector<uint8_t> drop;

  void set_constraint(size_t i)
  {
    const double l = A_pos[i] - A_pos[i - 1];
    A_b[i] = (A_cumsum[i] - A_cumsum[i - 1]) / l;
    A_q[i] = (A_cost[i] - A_cost[i - 1]) / l;
  }

  void push(double sum, double cost, double pos, bool constrained)
  {
    A_cumsum.push_back(sum);
    A_cost.push_back(cost);
    A_pos.push_back(pos);
    if (!constrained) return;
    A_b.push_back(0.0);
    A_q.push_back(0.0);
    if (A_pos.size() > 1) set_constraint(A_pos.size() - 1);
  }
};

/// VARIANCE model: + rounding errors of the sums (see 1D_Variance.h)
struct VarianceState : State
{
  std::vector<double> cumsum_lo, A_cumsum_lo;

  void set_constraint_variance(size_t i)
  {
    const double l = A_pos[i] - A_pos[i - 1];
    A_b[i] = ((A_cumsum[i] - A_cumsum[i - 1]) + (A_cumsum_lo[i] - A_cumsum_lo[i - 1])) / l;
    A_q[i] = (A_cost[i] - A_cost[i - 1]) / l;
  }

  void push_variance(double sum, double sum_lo, double cost, double pos, bool constrained)
  {
    A_cumsum.push_back(sum);
    A_cumsum_lo.push_back(sum_lo);
    A_cost.push_back(cost);
    A_pos.push_back(pos);
    if (!constrained) return;
    A_b.push_back(0.0);
    A_q.push_back(0.0);
    if (A_pos.size() > 1) set_constraint_variance(A_pos.size() - 1);
  }
};

} // namespace dust_hw

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "DUST_1D_HW-inl.h"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "DUST_1D_HW-inl.h"

namespace dust_hw {
HWY_EXPORT(Advance);
HWY_EXPORT(AdvanceVariance);
}

////////////////////////////////////////////////////////////////////////////////
/// Highway object (append_data, then update_partition from the last time)
class DUST_1D_HW_Obj
{
  public:
    DUST_1D_HW_Obj(std::string model, std::string method)
      : model_name(std::move(model)), method_name(std::move(method))
    {
      static const char* const methods[] = {"DUST", "PELT", "OP", "DUSTib"};
      for (int K = 0; K < 8; K++)
        dust_hw::with_model(K, [&](auto k) {
          if (model_name == dust_hw::ScalarModel<decltype(k)::value>::get_model()) s.model = K;
        });
      s.method = std::find(methods, methods + 4, method_name) - methods;
      if (s.model < 0) stop("Unrecognized model \"" + model_name + "\".");
      if (s.method == 4) stop("Unrecognized method \"" + method_name + "\". Valid values: \"DUST\", \"DUSTib\", \"PELT\", \"OP\".");
    }

    void append_data(const NumericVector& data, Nullable<double> inPenalty)
    {
      const size_t m = data.size();
      if (m == 0) return;
      dust_hw::with_model(s.model, [&](auto k) {
        for (double y : data) dustib::Math<decltype(k)::value>::validate(y);
      });
      const bool first_execution = s.n == 0;
      const double checked_penalty = dust_validation::penalty(inPenalty, s.penalty, m, !first_execution);
      s.n += m;

      if (first_execution)
      {
        s.penalty = checked_penalty;
        s.cumsum.push_back(0.0);
        s.costRecord.push_back(-s.penalty);
        s.chptRecord.push_back(0);
        if (variance()) { s.cumsum_lo.push_back(0.0); s.push_variance(0.0, 0.0, -s.penalty, 0.0, true); }
        else s.push(0.0, -s.penalty, 0.0, true);
      }
      s.cumsum.reserve(s.n + 1);
      s.costRecord.reserve(s.n + 1);
      s.chptRecord.reserve(s.n + 1);
      s.nb.reserve(s.n);
      if (variance())
      {
        // sums of y^2 + rounding errors (TwoSum)
        s.cumsum_lo.reserve(s.n + 1);
        for (double y : data)
        {
          const double x = y * y, a = s.cumsum.back(), sum = a + x, b = sum - a;
          s.cumsum.push_back(sum);
          s.cumsum_lo.push_back(s.cumsum_lo.back() + ((a - (sum - b)) + (x - b)));
        }
        return;
      }
      dust_hw::with_model(s.model, [&](auto k) {
        for (double y : data)
          s.cumsum.push_back(s.cumsum.back() + dust_hw::ScalarModel<decltype(k)::value>::statistic(y));
      });
    }

    void update_partition()
    {
      if (s.n == 0) return;
      if (variance()) HWY_DYNAMIC_DISPATCH(dust_hw::AdvanceVariance)(s);
      else HWY_DYNAMIC_DISPATCH(dust_hw::Advance)(s);
    }

    List get_partition()
    {
      if (s.n == 0) stop("append data before requesting a partition");
      if (s.chptRecord.size() != static_cast<size_t>(s.n) + 1)
        stop("update_partition before requesting a partition");
      std::vector<int> changepoints;
      for (int tau = s.chptRecord[s.n]; tau != 0; tau = s.chptRecord[tau]) changepoints.push_back(tau);
      std::reverse(changepoints.begin(), changepoints.end());
      changepoints.push_back(static_cast<int>(s.n));

      return List::create(
        _["changepoints"] = changepoints,
        _["lastIndexSet"] = std::vector<unsigned>(s.A_pos.rbegin(), s.A_pos.rend()),
        _["backend"] = "highway",
        _["nb"] = s.nb,
        _["costQ"] = std::vector<double>(s.costRecord.begin() + 1, s.costRecord.end())
      );
    }

    List get_info()
    {
      return List::create(
        _["backend"] = "highway",
        _["data_statistic"] = s.cumsum,
        _["data_length"] = s.n,
        _["current_penalty"] = s.penalty,
        _["model"] = model_name,
        _["pruning_algo"] = method_name
      );
    }

    List dust(const NumericVector& data, Nullable<double> inPenalty)
    {
      append_data(data, inPenalty);
      update_partition();
      return get_partition();
    }

  private:
    dust_hw::VarianceState s;
    std::string model_name, method_name;
    // variance: own kernels (DUST_1D_HW_Variance-inl.h)
    bool variance() const { return s.model == 7; }
};

#else // no Highway

DUST_1D* newModule1D(const std::string& model, const std::string& method); // _ModuleAssembly_1D.cpp

class DUST_1D_HW_Obj
{
  public:
    DUST_1D_HW_Obj(std::string, std::string)
    {
      stop("This dust build has no Highway backend. dust.object.1D.HW() "
           "should have redirected to dust.object.1D() automatically -- "
           "seeing this means that redirect was bypassed.");
    }
    void append_data(const NumericVector&, Nullable<double>) {}
    void update_partition() {}
    List get_partition() { return List(); }
    List get_info() { return List(); }
    List dust(const NumericVector&, Nullable<double>) { return List(); }
};

#endif // HAVE_HIGHWAY

//' @title Highway segmentation module
//'
//' @name DUST_1D_HW_Obj
//'
//' @description
//' Rcpp module for \code{\link{dust.object.1D}} with Highway.
//'
//' @export
RCPP_MODULE(DUSTHWMODULE1D)
{
  class_<DUST_1D_HW_Obj>("DUST_1D_HW_Obj")
    .constructor<std::string, std::string>()
    .method("append_data", &DUST_1D_HW_Obj::append_data)
    .method("update_partition", &DUST_1D_HW_Obj::update_partition)
    .method("get_partition", &DUST_1D_HW_Obj::get_partition)
    .method("get_info", &DUST_1D_HW_Obj::get_info)
    .method("dust", &DUST_1D_HW_Obj::dust)
  ;
}
//' DUST.1D.HW.backend
//'
//' @description Backend available in this installation: "highway" or "scalar".
//' @return "highway" or "scalar"
//' @keywords internal
// [[Rcpp::export(name = "DUST.1D.HW.backend")]]
CharacterVector DUST_1D_HW_backend()
{
#ifdef HAVE_HIGHWAY
  return CharacterVector::create("highway");
#else
  return CharacterVector::create("scalar");
#endif
}

//' DUST.1D.HW
//'
//' @description
//' Same as \code{\link{dust.1D}} with the Highway (SIMD) engine.
//' Use \code{dust.1D(..., backend = "highway")}.
//'
//' @param data a numeric vector
//' @param penalty penalty value, \code{2 log(n)} by default
//' @param model one of "gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance"
//' @param method one of "DUST", "DUSTib", "PELT", "OP"
//' @return a list with changepoints, lastIndexSet, backend, nb and costQ
//' @examples
//' y <- dataGenerator_1D(chpts = c(500, 1000), parameters = c(0, 1), type = "gauss")
//' dust.1D(y, backend = "highway")$changepoints
//' @keywords internal
// [[Rcpp::export(name = "DUST.1D.HW")]]
List DUST_1D_HW(NumericVector data, Nullable<double> penalty = R_NilValue,
                std::string model = "gauss", std::string method = "DUST")
{
  if (data.size() == 0) stop("data must be a nonempty numeric vector");
#ifdef HAVE_HIGHWAY
  return DUST_1D_HW_Obj(model, method).dust(data, penalty);
#else
  return std::unique_ptr<DUST_1D>(newModule1D(model, method))->dust(data, penalty);
#endif
}
