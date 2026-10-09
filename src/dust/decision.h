#ifndef DUST_DECISION_H
#define DUST_DECISION_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "common.h"
#include "costs.h"

namespace dust {

/// infinite outside the domain
template <Cost K>
double sum_conjugate(const std::vector<double>& z, double& scale)
{
  for (double value : z)
    if (!Family<K>::valid(value)) return std::numeric_limits<double>::infinity();
  double total = 0.0;
  scale = 0.0;
  for (double value : z)
  {
    const double term = Family<K>::conjugate(value);
    total += term;
    scale += std::abs(term);
  }
  return total;
}

constexpr double kInside = 1e-12;

/// for the searches only, D is evaluated at the true mean
template <Cost K>
double interior_mean(double z)
{
  if constexpr (K == Cost::GaussianMean) return z;
  else if constexpr (Family<K>::bounded) return std::min(1.0 - kInside, std::max(kInside, z));
  else return std::max(z, Family<K>::lower + kInside);
}

/// D(x) = - sum A*(a + Mx) - c - u.x, x >= 0
template <Cost K>
struct Decision
{
  size_t dimension;
  size_t constraints;
  std::vector<double> a;
  std::vector<double> matrix; // one column per constraint
  std::vector<double> u;
  double c;
  mutable std::vector<double> scratch, work_x, work_z, work_v;

  void resize(size_t d, size_t k)
  {
    dimension = d;
    constraints = k;
    a.resize(d);
    scratch.resize(d);
    matrix.resize(d * k);
    u.resize(k);
  }

  void mean(const std::vector<double>& x, std::vector<double>& z) const
  {
    z.assign(a.begin(), a.end());
    for (size_t j = 0; j < constraints; ++j)
      for (size_t i = 0; i < dimension; ++i) z[i] += matrix[j * dimension + i] * x[j];
  }

  double value(const std::vector<double>& x, double& scale) const
  {
    // x has to stay in the positive orthant
    for (double v : x)
      if (!std::isfinite(v) || v < 0.0)
      {
        scale = std::numeric_limits<double>::infinity();
        return -std::numeric_limits<double>::infinity();
      }
    mean(x, scratch);
    double linear = c;
    for (size_t j = 0; j < constraints; ++j) linear += u[j] * x[j];
    double conjugate_scale = 0.0;
    const double conjugate = sum_conjugate<K>(scratch, conjugate_scale);
    scale = std::abs(c) + conjugate_scale;
    for (size_t j = 0; j < constraints; ++j)
      scale += std::abs(u[j] * x[j]);
    return -conjugate - linear;
  }

  bool positive(const std::vector<double>& x) const
  {
    double scale = 0.0;
    const double score = value(x, scale);
    return dust::positive(score, scale);
  }

  // false outside the domain; trace of the curvature for the first step
  bool derivatives(const std::vector<double>& x,
                   std::vector<double>& gradient, double* trace = nullptr) const
  {
    for (size_t j = 0; j < constraints; ++j) gradient[j] = -u[j];
    if (trace) *trace = 0.0;
    for (size_t i = 0; i < dimension; ++i)
    {
      double z = a[i];
      double norm2 = 0.0;
      for (size_t j = 0; j < constraints; ++j)
      {
        const double b = matrix[j * dimension + i];
        z += b * x[j];
        norm2 += b * b;
      }
      if (!Family<K>::valid(z)) return false;
      if (norm2 == 0.0) continue;
      // curvature is infinite on the boundary
      if (trace && !Family<K>::boundary(z)) *trace += norm2 * Family<K>::curvature(z);
      z = interior_mean<K>(z);
      const double theta = Family<K>::theta(z);
      if (!std::isfinite(theta)) return false;
      for (size_t j = 0; j < constraints; ++j)
        gradient[j] -= matrix[j * dimension + i] * theta;
    }
    return std::all_of(gradient.begin(), gradient.end(),
                       [](double g) { return std::isfinite(g); });
  }

  void base_without(const std::vector<double>& x, size_t j,
                    std::vector<double>& base) const
  {
    std::copy(a.begin(), a.end(), base.begin());
    for (size_t k = 0; k < constraints; ++k)
      if (k != j)
        for (size_t i = 0; i < dimension; ++i)
          base[i] += matrix[k * dimension + i] * x[k];
  }

  double coordinate_value(const std::vector<double>& base, size_t j,
                          double at, double linear_base,
                          double scale_base, double& scale) const
  {
    for (size_t i = 0; i < dimension; ++i)
      scratch[i] = base[i] + matrix[j * dimension + i] * at;
    double conjugate_scale = 0.0;
    const double conjugate = sum_conjugate<K>(scratch, conjugate_scale);
    scale = scale_base + std::abs(u[j] * at) + conjugate_scale;
    return -conjugate - linear_base - u[j] * at;
  }

  std::pair<double, double> coordinate_slope_curvature(
    const std::vector<double>& base, size_t j, double at) const
  {
    double first = -u[j], second = 0.0;
    for (size_t i = 0; i < dimension; ++i)
    {
      const double m = matrix[j * dimension + i];
      if (m == 0.0) continue;
      double z = base[i] + m * at;
      z = interior_mean<K>(z);
      first -= m * Family<K>::theta(z);
      second -= m * m * Family<K>::curvature(z);
    }
    return {first, second};
  }

  std::pair<double, double> coordinate_bounds(const std::vector<double>& base,
                                               size_t j) const
  {
    const double lower = K == Cost::GaussianMean ? -std::numeric_limits<double>::infinity() : Family<K>::lower;
    const double upper = Family<K>::bounded ? 1.0 : std::numeric_limits<double>::infinity();
    double left = 0.0;
    double right = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < dimension; ++i)
    {
      const double slope = matrix[j * dimension + i];
      if (slope == 0.0) continue;
      if (slope > 0.0)
      {
        if (std::isfinite(lower)) left = std::max(left, (lower - base[i]) / slope);
        if (std::isfinite(upper)) right = std::min(right, (upper - base[i]) / slope);
      }
      else
      {
        if (std::isfinite(upper)) left = std::max(left, (upper - base[i]) / slope);
        if (std::isfinite(lower)) right = std::min(right, (lower - base[i]) / slope);
      }
    }
    if constexpr (K == Cost::Exponential || K == Cost::GaussianVariance)
    {
      if (left > 0.0) left = std::nextafter(left, std::numeric_limits<double>::infinity());
      if (std::isfinite(right)) right = std::nextafter(right, 0.0);
    }
    return {left, right};
  }
};

/// weak duality: max D <= A(v) - a.v - c for M'v + u >= 0, v = theta(z) projected
template <Cost K>
bool no_pruning_certificate(const Decision<K>& test, const double* z)
{
  const size_t d = test.dimension, k = test.constraints;
  std::vector<double>& v = test.work_v;
  v.resize(d);
  for (size_t i = 0; i < d; ++i)
  {
    if (!Family<K>::valid(z[i])) return false;
    v[i] = Family<K>::theta(interior_mean<K>(z[i]));
    if (!std::isfinite(v[i])) return false;
  }
  const auto slack = [&](size_t j) {
    double r = test.u[j];
    for (size_t i = 0; i < d; ++i) r += test.matrix[j * d + i] * v[i];
    return r;
  };
  // cyclic projections (one pass if one constraint)
  for (unsigned int pass = 0; pass < (k == 1 ? 1u : 8u); ++pass)
  {
    bool feasible = true;
    for (size_t j = 0; j < k; ++j)
    {
      const double r = slack(j);
      if (r >= 0.0) continue;
      feasible = false;
      double norm2 = 0.0;
      for (size_t i = 0; i < d; ++i)
        norm2 += test.matrix[j * d + i] * test.matrix[j * d + i];
      // D unbounded
      if (!(norm2 > 0.0)) return false;
      // a bit inside the half-space
      const double step = (kInside * (1.0 + std::abs(test.u[j])) - r) / norm2;
      for (size_t i = 0; i < d; ++i) v[i] += step * test.matrix[j * d + i];
    }
    if (feasible) break;
  }
  for (size_t j = 0; j < k; ++j)
    if (!(slack(j) >= 0.0)) return false;
  double bound = -test.c, scale = std::abs(test.c);
  for (size_t i = 0; i < d; ++i)
  {
    const double partition = Family<K>::partition(v[i]);
    if (!std::isfinite(partition)) return false; // v not in the domain
    bound += partition - test.a[i] * v[i];
    scale += std::abs(partition) + std::abs(test.a[i] * v[i]);
  }
  return std::isfinite(bound) && positive(-bound, scale);
}

/// argmax of D along column j on [lower, upper]: Newton, bisection or doubling; {point, unbounded}
template <Cost K>
std::pair<double, bool> coordinate_argmax(const Decision<K>& test,
                                          const std::vector<double>& base,
                                          size_t j, double lower, double upper)
{
  double lo = lower, hi = upper;
  double at = lower;
  auto derivatives = test.coordinate_slope_curvature(base, j, at);
  for (unsigned int k = 0; k < 64; ++k)
  {
    double next = at - derivatives.first / derivatives.second;
    if (!std::isfinite(next) || !(next > lo) || !(next < hi))
      next = std::isfinite(hi) ? lo + 0.5 * (hi - lo) :
        std::max(2.0 * at, lower + 1.0);
    if (!std::isfinite(next)) return {at, true};
    derivatives = test.coordinate_slope_curvature(base, j, next);
    if (derivatives.first > 0.0) lo = next; else hi = next;
    const bool converged = std::abs(next - at) <= 1e-10 * (1.0 + std::abs(next)) ||
      derivatives.first == 0.0;
    at = next;
    // converged (slope can stay > 0 when coming from the left)
    if (converged) return {at, false};
  }
  return {at, derivatives.first > 0.0 && !std::isfinite(hi)};
}

} // namespace dust

#endif
