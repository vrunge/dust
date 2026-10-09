/// Costs of a segment
///
/// The costs are exponential families with sufficient statistics per sample
/// and feature: y, or y^2 for GaussianVariance, (y, y^2) for
/// GaussianMeanVariance. For the one-parameter costs K and the mean a of the
/// statistic on a segment, the segment cost is -length * conjugate(a). The
/// pruning tests use the conjugate A*, its derivative theta = A*' (natural
/// parameter), mean = theta^-1, the log-partition A and the curvature A*''.

#ifndef DUST_COSTS_H
#define DUST_COSTS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace dust {

/// cost of a segment (Cost::GaussianMean: GaussianMeanCost in the Python package changepoints)
enum class Cost
{
  GaussianMean,       // y ~ N(mu, 1)
  Poisson,
  Exponential,
  Geometric,          // number of trials (>= 1)
  Bernoulli,          // y in [0, 1]
  Binomial,           // proportion of successes (y / m)
  NegativeBinomial,   // y / r
  GaussianVariance,   // y ~ N(0, sigma^2)
  GaussianMeanVariance
};

/// f(std::integral_constant<Cost, K>{}) for the cost K
template <class F>
decltype(auto) with_cost(Cost K, F&& f)
{
  switch (K)
  {
    case Cost::GaussianMean:     return f(std::integral_constant<Cost, Cost::GaussianMean>{});
    case Cost::Poisson:          return f(std::integral_constant<Cost, Cost::Poisson>{});
    case Cost::Exponential:      return f(std::integral_constant<Cost, Cost::Exponential>{});
    case Cost::Geometric:        return f(std::integral_constant<Cost, Cost::Geometric>{});
    case Cost::Bernoulli:        return f(std::integral_constant<Cost, Cost::Bernoulli>{});
    case Cost::Binomial:         return f(std::integral_constant<Cost, Cost::Binomial>{});
    case Cost::NegativeBinomial: return f(std::integral_constant<Cost, Cost::NegativeBinomial>{});
    case Cost::GaussianVariance: return f(std::integral_constant<Cost, Cost::GaussianVariance>{});
    default:                     return f(std::integral_constant<Cost, Cost::GaussianMeanVariance>{});
  }
}

template <Cost K>
struct Family
{
  static constexpr size_t statistics = K == Cost::GaussianMeanVariance ? 2 : 1;   // per feature
  static constexpr bool bounded = K == Cost::Bernoulli || K == Cost::Binomial;    // mean in [0, 1]
  static constexpr bool negative = K == Cost::Exponential || K == Cost::Geometric ||
                                   K == Cost::NegativeBinomial || K == Cost::GaussianVariance;   // theta < 0
  static constexpr double lower = K == Cost::Geometric ? 1.0 : 0.0;   // smallest mean (not gauss)

  static void statistic(double y, double* s)
  {
    s[0] = K == Cost::GaussianVariance ? y * y : y;
    if constexpr (statistics == 2) s[1] = y * y;
  }

  static bool valid(double a)
  {
    if (!std::isfinite(a)) return false;
    if constexpr (K == Cost::GaussianMean) return true;
    if constexpr (bounded) return a >= 0 && a <= 1;
    if constexpr (K == Cost::Exponential || K == Cost::GaussianVariance) return a > 0;
    return a >= lower;
  }

  /// mean on a boundary of the domain (one-point segment)
  static bool boundary(double a)
  {
    if constexpr (K == Cost::GaussianMean || K == Cost::Exponential || K == Cost::GaussianVariance) return false;
    if constexpr (bounded) return a == 0 || a == 1;
    return a == lower;
  }

  static void validate(double y)
  {
    if (K == Cost::GaussianMeanVariance ? std::isfinite(y) : valid(K == Cost::GaussianVariance ? y * y : y)) return;
    throw std::domain_error("observation outside the model domain (finite values; Exponential and Variance "
                            "require strictly positive sufficient statistics; Binomial requires proportions in [0,1])");
  }

  static double conjugate(double a)
  {
    if (boundary(a)) return 0;
    if constexpr (K == Cost::GaussianMean) return .5 * a * a;
    if constexpr (K == Cost::Poisson) return a * (std::log(a) - 1);
    if constexpr (K == Cost::Exponential) return -std::log(a) - 1;
    if constexpr (K == Cost::GaussianVariance) return -.5 * (std::log(a) + 1);
    if constexpr (bounded) return a * std::log(a) + (1 - a) * std::log1p(-a);
    if constexpr (K == Cost::Geometric) return (a - 1) * std::log1p(-1 / a) - std::log(a);
    if constexpr (K == Cost::NegativeBinomial)
    {
      // (log1p(-1/(1+a)) not accurate for small a)
      const double l = a < 1 ? std::log(a) - std::log1p(a) : -std::log1p(1 / a);
      return a * l - std::log1p(a);
    }
  }

  static double theta(double a)
  {
    if constexpr (K == Cost::GaussianMean) return a;
    if constexpr (K == Cost::Poisson) return std::log(a);
    if constexpr (K == Cost::Exponential) return -1 / a;
    if constexpr (K == Cost::GaussianVariance) return -.5 / a;
    if constexpr (bounded) return std::log(a) - std::log1p(-a);
    if constexpr (K == Cost::Geometric) return std::log1p(-1 / a);
    if constexpr (K == Cost::NegativeBinomial) return a < 1 ? std::log(a) - std::log1p(a) : -std::log1p(1 / a);
  }

  static double mean(double r)
  {
    if constexpr (K == Cost::GaussianMean) return r;
    if constexpr (K == Cost::Poisson) return std::exp(r);
    if constexpr (K == Cost::Exponential) return -1 / r;
    if constexpr (K == Cost::GaussianVariance) return -.5 / r;
    if constexpr (bounded) { const double z = std::exp(-std::abs(r)); return r >= 0 ? 1 / (1 + z) : z / (1 + z); }
    if constexpr (K == Cost::Geometric) return -1 / std::expm1(r);
    if constexpr (K == Cost::NegativeBinomial) return std::exp(r) / -std::expm1(r);
  }

  static double partition(double r)
  {
    if constexpr (K == Cost::GaussianMean) return .5 * r * r;
    if constexpr (K == Cost::Poisson) return std::exp(r);
    if constexpr (K == Cost::Exponential) return -std::log(-r);
    if constexpr (K == Cost::GaussianVariance) return -.5 * (std::log(-r) + std::log(2.0));
    if constexpr (bounded) return std::max(r, 0.0) + std::log1p(std::exp(-std::abs(r)));
    if constexpr (K == Cost::Geometric || K == Cost::NegativeBinomial)
    {
      const double l = r < -std::log(2.0) ? std::log1p(-std::exp(r)) : std::log(-std::expm1(r));
      return (K == Cost::Geometric ? r : 0) - l;
    }
  }

  /// A*''
  static double curvature(double a)
  {
    if constexpr (K == Cost::GaussianMean) return 1;
    if constexpr (K == Cost::Poisson) return 1 / a;
    if constexpr (K == Cost::Exponential) return 1 / (a * a);
    if constexpr (K == Cost::GaussianVariance) return .5 / (a * a);
    if constexpr (bounded) return 1 / a + 1 / (1 - a);
    if constexpr (K == Cost::Geometric) return 1 / (a - 1) - 1 / a;
    if constexpr (K == Cost::NegativeBinomial) return 1 / (a * (1 + a));
  }
};

} // namespace dust

#endif
