/// MD DUST: the decision function (included by MD_DUST.cpp)

#ifndef MD_DECISION_H
#define MD_DECISION_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "1D_DUSTib.h"
#ifdef HAVE_HIGHWAY
#include "DUST_1D_HW_Models.h"
#endif

namespace dust_md {

using namespace Rcpp;

template <class Model> using Math = dustib::Math<dustib::model_id<Model>>;

#ifdef HAVE_HIGHWAY
template <class Model> struct VectorModel;
template <> struct VectorModel<GaussPolicy>    { using type = VecGaussPolicy; };
template <> struct VectorModel<PoissonPolicy>  { using type = VecPoissonPolicy; };
template <> struct VectorModel<ExpPolicy>      { using type = VecExpPolicy; };
template <> struct VectorModel<GeomPolicy>     { using type = VecGeomPolicy; };
template <> struct VectorModel<BernPolicy>     { using type = VecBernPolicy; };
template <> struct VectorModel<BinomPolicy>    { using type = VecBinomPolicy; };
template <> struct VectorModel<NegbinPolicy>   { using type = VecNegbinPolicy; };
template <> struct VectorModel<VariancePolicy> { using type = VecVariancePolicy; };
#endif

template <class Model>
double sum_conjugate(const std::vector<double>& z, bool highway, double& scale)
{
  constexpr int K = dustib::model_id<Model>;
  for (double value : z)
    if (!Math<Model>::valid(value)) return std::numeric_limits<double>::infinity();

  double total = 0.0;
  scale = 0.0;
#ifdef HAVE_HIGHWAY
  if (highway && z.size() >= hn::Lanes(hn::ScalableTag<double>()))
  {
    bool interior = true;
    for (double value : z)
    {
      if (Math<Model>::boundary(value) ||
          ((K == 1 || K == 6 || K == 2 || K == 7) && value < 1e-8) ||
          (K == 3 && value < 1.0 + 1e-8) ||
          ((K == 4 || K == 5) && (value < 1e-8 || value > 1.0 - 1e-8)))
      {
        interior = false;
        break;
      }
    }
    if (interior)
    {
      using Vec = typename VectorModel<Model>::type;
      const hn::ScalableTag<double> tag;
      const size_t lanes = hn::Lanes(tag);
      double values[hn::MaxLanes(tag)];
      size_t i = 0;
      for (; i + lanes <= z.size(); i += lanes)
      {
        hn::StoreU(Vec::vDstar(tag, hn::LoadU(tag, z.data() + i)), tag, values);
        for (size_t j = 0; j < lanes; ++j)
        {
          total += values[j];
          scale += std::abs(values[j]);
        }
      }
      for (; i < z.size(); ++i)
      {
        const double value = Math<Model>::conjugate(z[i]);
        total += value;
        scale += std::abs(value);
      }
      return total;
    }
  }
#else
  (void)highway;
  (void)K;
#endif
  for (double value : z)
  {
    const double term = Math<Model>::conjugate(value);
    total += term;
    scale += std::abs(term);
  }
  return total;
}

/// mean moved inside the domain (for the searches only, D is always
/// evaluated at the true mean)
template <class Model>
double interior_mean(double z)
{
  if constexpr (dustib::model_id<Model> == 3)
    return std::max(z, 1.0 + 1e-12);
  else if constexpr (dustib::model_id<Model> == 4 || dustib::model_id<Model> == 5)
    return std::min(1.0 - 1e-12, std::max(1e-12, z));
  else if constexpr (dustib::model_id<Model> != 0)
    return std::max(z, 1e-12);
  else
    return z;
}

////////////////////////////////////////////////////////////////////////////////
/// decision function D(x) = - sum A*(a + Mx) - c - u.x, x >= 0
/// (Runge, Truong and Querne 2025, arXiv:2507.02467)
template <class Model>
struct Decision
{
  size_t dimension;
  size_t constraints;
  bool highway;
  std::vector<double> a;
  std::vector<double> matrix; // one column per constraint
  std::vector<double> u;
  double c;
  mutable std::vector<double> scratch;
  mutable std::vector<double> work_x, work_z, work_norm; // workspace
  mutable std::vector<double> work_v;

  double value(const std::vector<double>& x, double& scale) const
  {
    // x has to stay in the positive orthant
    for (double v : x)
      if (!std::isfinite(v) || v < 0.0)
      {
        scale = std::numeric_limits<double>::infinity();
        return -std::numeric_limits<double>::infinity();
      }
    std::copy(a.begin(), a.end(), scratch.begin());
    double linear = c;
    for (size_t j = 0; j < constraints; ++j)
    {
      linear += u[j] * x[j];
      for (size_t i = 0; i < dimension; ++i)
        scratch[i] += matrix[j * dimension + i] * x[j];
    }
    double conjugate_scale = 0.0;
    const double conjugate = sum_conjugate<Model>(scratch, highway, conjugate_scale);
    scale = std::abs(c) + conjugate_scale;
    for (size_t j = 0; j < constraints; ++j)
      scale += std::abs(u[j] * x[j]);
    return -conjugate - linear;
  }

  bool positive(const std::vector<double>& x) const
  {
    double scale = 0.0;
    const double score = value(x, scale);
    return dustib::positive(score, scale);
  }

  // gradient of D at x (+ trace of the curvature for the first step)
  // false if x is outside the domain
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
      if (!Math<Model>::valid(z)) return false;
      if (norm2 == 0.0) continue;
      // curvature is infinite on the boundary
      if (trace && !Math<Model>::boundary(z)) *trace += norm2 * Model::DstarSecond(z);
      z = interior_mean<Model>(z);
      const double theta = Math<Model>::theta(z);
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
    const double conjugate = sum_conjugate<Model>(scratch, highway, conjugate_scale);
    scale = scale_base + std::abs(u[j] * at) + conjugate_scale;
    return -conjugate - linear_base - u[j] * at;
  }

  bool coordinate_positive(const std::vector<double>& base, size_t j,
                           double at, double linear_base,
                           double scale_base) const
  {
    double scale = 0.0;
    return dustib::positive(
      coordinate_value(base, j, at, linear_base, scale_base, scale), scale);
  }

  double coordinate_derivative(const std::vector<double>& base,
                               size_t j, double at) const
  {
    double result = -u[j];
    for (size_t i = 0; i < dimension; ++i)
    {
      double z = base[i] + matrix[j * dimension + i] * at;
      z = interior_mean<Model>(z);
      result -= matrix[j * dimension + i] * Math<Model>::theta(z);
    }
    return result;
  }

  // first and second derivative of D along coordinate j
  std::pair<double, double> coordinate_slope_curvature(
    const std::vector<double>& base, size_t j, double at) const
  {
    double first = -u[j], second = 0.0;
    for (size_t i = 0; i < dimension; ++i)
    {
      const double m = matrix[j * dimension + i];
      if (m == 0.0) continue;
      double z = base[i] + m * at;
      z = interior_mean<Model>(z);
      first -= m * Math<Model>::theta(z);
      second -= m * m * Model::DstarSecond(z);
    }
    return {first, second};
  }

  std::pair<double, double> coordinate_bounds(const std::vector<double>& base,
                                               size_t j) const
  {
    constexpr int K = dustib::model_id<Model>;
    const double lower = K == 3 ? 1.0 :
      K == 0 ? -std::numeric_limits<double>::infinity() : 0.0;
    const double upper = K == 4 || K == 5 ? 1.0 :
      std::numeric_limits<double>::infinity();
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
    if constexpr (K == 2 || K == 7)
    {
      if (left > 0.0) left = std::nextafter(left, std::numeric_limits<double>::infinity());
      if (std::isfinite(right)) right = std::nextafter(right, 0.0);
    }
    return {left, right};
  }
};


////////////////////////////////////////////////////////////////////////////////
/// NO PRUNING certificate (weak duality)
/// for v with M'v + u >= 0 : max D <= A(v) - a.v - c
/// v = theta(z) at the current point, projected on {M'v + u >= 0}
/// true if this bound is negative -> we can stop the search
template <class Model>
bool no_pruning_certificate(const Decision<Model>& test, const double* z)
{
  const size_t d = test.dimension, k = test.constraints;
  std::vector<double>& v = test.work_v;
  v.resize(d);
  for (size_t i = 0; i < d; ++i)
  {
    if (!Math<Model>::valid(z[i])) return false;
    v[i] = Math<Model>::theta(interior_mean<Model>(z[i]));
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
      const double step = (1e-12 * (1.0 + std::abs(test.u[j])) - r) / norm2;
      for (size_t i = 0; i < d; ++i) v[i] += step * test.matrix[j * d + i];
    }
    if (feasible) break;
  }
  for (size_t j = 0; j < k; ++j)
    if (!(slack(j) >= 0.0)) return false;
  double bound = -test.c, scale = std::abs(test.c);
  for (size_t i = 0; i < d; ++i)
  {
    const double partition = Math<Model>::partition(v[i]);
    if (!std::isfinite(partition)) return false; // v not in the domain
    bound += partition - test.a[i] * v[i];
    scale += std::abs(partition) + std::abs(test.a[i] * v[i]);
  }
  return std::isfinite(bound) && dustib::positive(-bound, scale);
}

////////////////////////////////////////////////////////////////////////////////
/// argmax of D along column j on [lower, upper] (positive slope at lower)
/// Newton steps inside the bracket, bisection or doubling otherwise
/// returns {point, unbounded}
template <class Model>
std::pair<double, bool> coordinate_argmax(const Decision<Model>& test,
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

} // namespace dust_md

#endif
