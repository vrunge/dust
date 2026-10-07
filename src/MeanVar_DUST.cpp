#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();

/// Gaussian cost (mean and variance), infinite for 1 point or variance 0
double segment_cost(double sum, double sum2, double len)
{
  if (len < 2.0) return inf;
  const double mean = sum / len;
  const double variance = sum2 / len - mean * mean;
  if (!(variance > 0.0) || !std::isfinite(variance)) return inf;
  return 0.5 * len * (1.0 + std::log(variance));
}

} // namespace

#ifdef HAVE_HIGHWAY
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "MeanVar_DUST-inl.h"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "MeanVar_DUST-inl.h"
namespace meanvar_hw {
HWY_EXPORT(Scan);
}
#endif

namespace {

struct DecisionData {
  double m, m2, q;
};

/// one constraint: max of D(x, 0) (with x = 0)
template<bool check_zero_first = false>
bool one_constraint(const DecisionData& a, const DecisionData& b)
{
  if (!std::isfinite(a.q) || !std::isfinite(b.q)) {
    // infinite cost: test of the original meanVar code
    const double va = a.m2 - std::pow(a.m, 2);
    const double vb = b.m2 - std::pow(b.m, 2);
    const double delta2 = std::pow(a.m - b.m, 2);
    const double x0 = 0.5 * ((va - vb) / delta2 - 1.0);
    const double x1 = x0 * x0 + va / delta2;
    const double x2 = a.q - b.q;
    const double sign = x2 > 0.0 ? 1.0 : (x2 < 0.0 ? -1.0 : 0.0);
    const double x = std::max(0.0, x0 + 1.0 / (2.0 * x2)
      - sign * std::sqrt(x1 + 1.0 / (4.0 * x2 * x2)));
    const double A = a.m2 + x * (a.m2 - b.m2);
    const double B = a.m + x * (a.m - b.m);
    return 0.5 * (1.0 + std::log(A - B * B))
      - (a.q + x * (a.q - b.q)) > 0.0;
  }
  const double va = a.m2 - a.m * a.m;
  const double vb = b.m2 - b.m * b.m;
  if (!(va > 0.0) || !std::isfinite(va) || !std::isfinite(vb)) return false;
  const double delta = a.m - b.m;
  const double delta2 = delta * delta;
  const double linear = va - vb - delta2;
  const double slope = a.q - b.q;
  auto evaluate = [&](double x) {
    if (!(x >= 0.0) || !std::isfinite(x)) return false;
    const double v = va + x * linear - x * x * delta2;
    if (!(v > 0.0) || !std::isfinite(v)) return false;
    const double value = 0.5 * (1.0 + std::log(v)) - a.q - x * slope;
    return std::isfinite(value) && value > 0.0;
  };
  if constexpr (check_zero_first) {
    if (evaluate(0.0)) return true;
  }
  if (delta2 == 0.0) {
    // same means: variance linear in x
    if (linear >= 0.0 && slope < 0.0) return true;
    if (linear > 0.0 && slope == 0.0) return true;
    if (linear == 0.0 || slope == 0.0) return evaluate(0.0);
    return evaluate(std::max(0.0, 1.0 / (2.0 * slope) - va / linear));
  }
  const double x0 = linear / (2.0 * delta2);
  if (slope == 0.0) return evaluate(std::max(0.0, x0));

  // critical point (dust::decisionTest1)
  const double root = x0 * x0 + va / delta2 + 1.0 / (4.0 * slope * slope);
  if (!(root >= 0.0) || !std::isfinite(root)) return evaluate(0.0);
  const double sign = slope > 0.0 ? 1.0 : -1.0;
  const double x = x0 + 1.0 / (2.0 * slope) - sign * std::sqrt(root);
  return evaluate(std::max(0.0, x));
}

double decision_value(const DecisionData& a, const DecisionData& b1,
                      const DecisionData& b2, double x1, double x2)
{
  if (!(x1 >= 0.0 && x2 >= 0.0) || !std::isfinite(x1) || !std::isfinite(x2)) return -inf;
  const double mean = a.m + x1 * (a.m - b1.m) + x2 * (a.m - b2.m);
  const double second = a.m2 + x1 * (a.m2 - b1.m2) + x2 * (a.m2 - b2.m2);
  const double variance = second - mean * mean;
  if (!(variance > 0.0) || !std::isfinite(variance)) return -inf;
  return 0.5 * (1.0 + std::log(variance)) - a.q
    - x1 * (a.q - b1.q) - x2 * (a.q - b2.q);
}

bool two_constraints(const DecisionData& a, const DecisionData& b1,
                     const DecisionData& b2)
{
  if (one_constraint(a, b1) || one_constraint(a, b2)) return true;
  if (!std::isfinite(a.q) || !std::isfinite(b1.q) || !std::isfinite(b2.q)) return false;

  const double a0 = a.m2, a1 = a.m2 - b1.m2, a2 = a.m2 - b2.m2;
  const double b0 = a.m, b_1 = a.m - b1.m, b_2 = a.m - b2.m;
  const double c1 = a.q - b1.q, c2 = a.q - b2.q;
  const double denomA = 2.0 * (b_1 * c2 - b_2 * c1);
  if (!std::isfinite(denomA) || std::abs(denomA) < 1e-12)
    return decision_value(a, b1, b2, 0.0, 0.0) > 0.0;
  const double A = (b_1 * a2 - b_2 * a1) / denomA;
  double y;
  if (std::abs(b_1) >= std::abs(b_2) && std::abs(b_1) > 1e-12)
    y = (a1 - 2.0 * c1 * A) / (2.0 * b_1);
  else if (std::abs(b_2) > 1e-12)
    y = (a2 - 2.0 * c2 * A) / (2.0 * b_2);
  else return decision_value(a, b1, b2, 0.0, 0.0) > 0.0;

  const double k1 = A + y * y - a0;
  const double k2 = y - b0;
  const double det = a1 * b_2 - a2 * b_1;
  if (!std::isfinite(det) || std::abs(det) < 1e-12)
    return decision_value(a, b1, b2, 0.0, 0.0) > 0.0;
  const double x1 = (k1 * b_2 - k2 * a2) / det;
  const double x2 = (a1 * k2 - b_1 * k1) / det;
  return decision_value(a, b1, b2, x1, x2) > 0.0;
}

class MeanVarDust {
public:
  MeanVarDust(std::string method, std::string requested_backend)
    : method_(std::move(method))
  {
    if (method_ != "1D" && method_ != "2D" && method_ != "PELT")
      stop("method must be '1D', '2D', or 'PELT'");
    if (requested_backend != "highway" && requested_backend != "scalar")
      stop("backend must be 'highway' or 'scalar'");
#ifdef HAVE_HIGHWAY
    highway_ = requested_backend == "highway";
#endif
  }

  void append_data(NumericVector data, Nullable<double> inPenalty)
  {
    for (double y : data)
      if (!std::isfinite(y)) stop("meanVar data must contain only finite numbers");
    if (!inPenalty.isNull()) {
      const double supplied = as<double>(inPenalty);
      if (!std::isfinite(supplied) || supplied < 0.0)
        stop("penalty must be finite and nonnegative");
      if (initialized_ && supplied != penalty_)
        stop("penalty cannot change after the first nonempty append");
    }
    if (data.size() == 0) return;
    if (!initialized_) {
      penalty_ = inPenalty.isNull() ? 4.0 * std::log(static_cast<double>(data.size())) : as<double>(inPenalty);
      sums_.push_back(0.0);
      sums2_.push_back(0.0);
      sums_lo_.push_back(0.0);
      sums2_lo_.push_back(0.0);
      costs_.push_back(-penalty_);
      changepoints_.push_back(0);
      active_.push_back(0);
      if (highway_) {
        active_sums_.push_back(0.0);
        active_sums2_.push_back(0.0);
        active_sums_lo_.push_back(0.0);
        active_sums2_lo_.push_back(0.0);
        active_costs_.push_back(-penalty_);
        active_positions_.push_back(0.0);
      }
      initialized_ = true;
    }
    // sums of y and y^2 + rounding errors (see 1D_Variance.h)
    for (double y : data) {
      add(sums_, sums_lo_, y);
      add(sums2_, sums2_lo_, y * y);
    }
  }

  void update_partition()
  {
    if (!initialized_) return;
    for (uint32_t t = static_cast<uint32_t>(costs_.size()); t < sums_.size(); ++t) {
      const auto best = highway_ ? scan_highway(t) : scan_scalar(t);
      const double qt = best.first + penalty_;
      costs_.push_back(qt);
      changepoints_.push_back(best.second);

      const uint32_t previous_t = t - 1;
      const double previous_q = costs_[previous_t];
      drop_.assign(active_.size(), 0);
      bool any_pruned = false;
      for (size_t i = active_.size(); i-- > 1;) {
        const uint32_t s = active_[i], r1 = active_[i - 1];
        if (previous_t <= s + 1) continue;
        const uint32_t r2 = i >= 2 ? active_[i - 2] : r1;
        const DecisionData a = stats(s, previous_t, previous_q);
        const DecisionData b1 = stats(r1, s, costs_[s]);
        bool prune;
        if (method_ == "PELT") {
          const double cost = segment_cost(sum1(previous_t, s), sum2(previous_t, s),
                                           previous_t - s);
          // no pruning with an infinite cost (segment can become valid)
          prune = std::isfinite(cost) && previous_q < inf &&
            costs_[s] + cost > previous_q;
        } else {
          prune = method_ == "1D" ? one_constraint<true>(a, b1)
            : two_constraints(a, b1, stats(r2, s, costs_[s]));
        }
        drop_[i] = static_cast<uint8_t>(prune);
        any_pruned |= prune;
      }

      if (any_pruned) {
        size_t write = 0;
        for (size_t i = 0; i < active_.size(); ++i) {
          if (drop_[i]) continue;
          active_[write] = active_[i];
          if (highway_) {
            active_sums_[write] = active_sums_[i];
            active_sums2_[write] = active_sums2_[i];
            active_sums_lo_[write] = active_sums_lo_[i];
            active_sums2_lo_[write] = active_sums2_lo_[i];
            active_costs_[write] = active_costs_[i];
            active_positions_[write] = active_positions_[i];
          }
          ++write;
        }
        active_.resize(write);
        if (highway_) {
          active_sums_.resize(write);
          active_sums2_.resize(write);
          active_sums_lo_.resize(write);
          active_sums2_lo_.resize(write);
          active_costs_.resize(write);
          active_positions_.resize(write);
        }
      }
      counts_.push_back(static_cast<int>(active_.size()));
      active_.push_back(t);
      if (highway_) {
        active_sums_.push_back(sums_[t]);
        active_sums2_.push_back(sums2_[t]);
        active_sums_lo_.push_back(sums_lo_[t]);
        active_sums2_lo_.push_back(sums2_lo_[t]);
        active_costs_.push_back(qt);
        active_positions_.push_back(static_cast<double>(t));
      }
    }
  }

  List get_partition() const
  {
    if (!initialized_) stop("append data before requesting a partition");
    if (costs_.size() != sums_.size())
      stop("update_partition before requesting a partition");
    const uint32_t n = static_cast<uint32_t>(costs_.size() - 1);
    if (n == 0 || !std::isfinite(costs_[n]))
      stop("no finite mean-variance segmentation is available for these data");
    std::vector<uint32_t> cp;
    for (uint32_t s = changepoints_[n]; s != 0; s = changepoints_[s]) cp.push_back(s);
    std::reverse(cp.begin(), cp.end());
    cp.push_back(n);
    return List::create(
      _["changepoints"] = cp,
      _["lastIndexSet"] = std::vector<uint32_t>(active_.rbegin(), active_.rend()),
      _["backend"] = highway_ ? "highway" : "scalar",
      _["nb"] = counts_,
      _["costQ"] = NumericVector(costs_.begin() + 1, costs_.end()));
  }

  List get_info() const
  {
    return List::create(
      _["backend"] = highway_ ? "highway" : "scalar",
      _["data_length"] = sums_.empty() ? 0 : sums_.size() - 1,
      _["current_penalty"] = penalty_,
      _["pruning_algo"] = method_,
      _["data_statistic"] = combined(sums_, sums_lo_),
      _["data_statistic2"] = combined(sums2_, sums2_lo_));
  }

  List dust(NumericVector data, Nullable<double> penalty)
  {
    append_data(data, penalty);
    update_partition();
    return get_partition();
  }

private:
  // TwoSum
  static void add(std::vector<double>& sums, std::vector<double>& lo, double x)
  {
    const double a = sums.back(), sum = a + x, b = sum - a;
    sums.push_back(sum);
    lo.push_back(lo.back() + ((a - (sum - b)) + (x - b)));
  }

  static std::vector<double> combined(const std::vector<double>& sums,
                                      const std::vector<double>& lo)
  {
    std::vector<double> out(sums.size());
    for (size_t i = 0; i < sums.size(); ++i) out[i] = sums[i] + lo[i];
    return out;
  }

  // sums of y and y^2 over (s, t]
  double sum1(uint32_t t, uint32_t s) const { return (sums_[t] - sums_[s]) + (sums_lo_[t] - sums_lo_[s]); }
  double sum2(uint32_t t, uint32_t s) const { return (sums2_[t] - sums2_[s]) + (sums2_lo_[t] - sums2_lo_[s]); }

  DecisionData stats(uint32_t s, uint32_t t, double qt) const
  {
    const double len = static_cast<double>(t - s);
    return {sum1(t, s) / len,
            sum2(t, s) / len,
            (qt - costs_[s]) / len};
  }

  std::pair<double, uint32_t> scan_scalar(uint32_t t) const
  {
    double best = inf;
    uint32_t index = 0;
    for (size_t i = active_.size(); i-- > 0;) {
      const uint32_t s = active_[i];
      const double value = costs_[s] + segment_cost(sum1(t, s), sum2(t, s), t - s);
      if (value < best) { best = value; index = s; }
    }
    return {best, index};
  }

  std::pair<double, uint32_t> scan_highway(uint32_t t) const
  {
#ifdef HAVE_HIGHWAY
    return HWY_DYNAMIC_DISPATCH(meanvar_hw::Scan)(
      active_positions_.data(), active_sums_.data(), active_sums_lo_.data(),
      active_sums2_.data(), active_sums2_lo_.data(), active_costs_.data(),
      active_.data(), active_.size(), t, sums_[t], sums_lo_[t], sums2_[t], sums2_lo_[t]);
#else
    return scan_scalar(t);
#endif
  }

  std::string method_;
  bool highway_ = false, initialized_ = false;
  double penalty_ = 0.0;
  std::vector<double> sums_, sums2_, costs_;
  std::vector<double> sums_lo_, sums2_lo_;   // rounding errors
  std::vector<double> active_sums_, active_sums2_, active_costs_, active_positions_;
  std::vector<double> active_sums_lo_, active_sums2_lo_;
  std::vector<uint8_t> drop_;
  std::vector<uint32_t> active_, changepoints_;
  std::vector<int> counts_;
};

MeanVarDust* new_meanvar(const std::string& method, const std::string& backend)
{
  return new MeanVarDust(method, backend);
}

} // namespace

RCPP_MODULE(DUSTMODULEmeanVar2)
{
  class_<MeanVarDust>("DUST_meanVar2")
    .factory<const std::string&, const std::string&>(new_meanvar)
    .method("append_data", &MeanVarDust::append_data)
    .method("update_partition", &MeanVarDust::update_partition)
    .method("get_partition", &MeanVarDust::get_partition)
    .method("get_info", &MeanVarDust::get_info)
    .method("dust", &MeanVarDust::dust);
}
