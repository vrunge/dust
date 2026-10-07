/// MD DUST: method "exact", D evaluated at its maximum
///
/// 1 constraint : one evaluation (closed form for gauss and for d = 1,
///                Newton steps otherwise)
/// 2 constraints: as in dust.meanVar, the 2 one-constraint maxima and the
///                critical point (exact for gauss, d = 2 otherwise)
/// k >= 3       : face solver for gauss (MD_GaussExact.h), 2 constraints
///                for the other models

#ifndef MD_EXACT_H
#define MD_EXACT_H

#include "MD_Decision.h"
#include "MD_GaussExact.h"

namespace dust_md {

////////////////////////////////////////////////////////////////////////////////
/// argmax of D in x_j (other multipliers = 0), returns {x_j, unbounded}
template <class Model>
std::pair<double, bool> single_constraint_argmax(const Decision<Model>& test,
                                                 size_t j)
{
  const size_t d = test.dimension;
  const double* m = test.matrix.data() + j * d;
  if constexpr (dustib::model_id<Model> == 0)
  {
    // D(x) = D(0) + g x - |m|^2 x^2 / 2, g = -u - m.a
    double norm = 0.0, g = -test.u[j];
    for (size_t i = 0; i < d; ++i)
    {
      norm += m[i] * m[i];
      g -= m[i] * test.a[i];
    }
    if (norm == 0.0) return {g > 0.0 ? 1.0 : 0.0, g > 0.0};
    return {std::max(0.0, g / norm), false};
  }
  else
  {
    const auto bounds = test.coordinate_bounds(test.a, j);
    const double lower = bounds.first, upper = bounds.second;
    if (!(upper > lower) || !(test.coordinate_derivative(test.a, j, lower) > 0.0))
      return {lower, false};
    if (d == 1 && m[0] != 0.0)
    {
      // closed form: theta(a + x m) = -u / m
      const double z = Math<Model>::mean(-test.u[j] / m[0]);
      const double x = (z - test.a[0]) / m[0];
      if (Math<Model>::valid(z) && std::isfinite(x))
        return {std::min(upper, std::max(lower, x)), false};
    }
    return coordinate_argmax(test, test.a, j, lower, upper);
  }
}

template <class Model>
bool positive_on_axis(const Decision<Model>& test, size_t j, double value)
{
  std::vector<double>& x = test.work_x;
  x.assign(test.constraints, 0.0);
  x[j] = value;
  return test.positive(x);
}

/// one constraint: evaluation at the maximum
template <class Model>
bool exact_one_constraint(const Decision<Model>& test, size_t j)
{
  const auto best = single_constraint_argmax(test, j);
  if (best.second)
  {
    // D unbounded: we go further until D > 0
    for (double value = 1.0; value < 1e300; value *= 2.0)
      if (positive_on_axis(test, j, value)) return true;
    return false;
  }
  return std::isfinite(best.first) && positive_on_axis(test, j, best.first);
}

////////////////////////////////////////////////////////////////////////////////
/// 2 constraints, D without maximum (not gauss): D can go to +infinity along
/// x = l (1 - w, w). We need M d >= 0 (M d = 0 for poisson, bern and binom),
/// then D ~ -u.d l and we try the two ends of the w interval.
template <class Model>
bool exact_two_unbounded(const Decision<Model>& test)
{
  constexpr int K = dustib::model_id<Model>;
  constexpr bool bounded = dustib::Math<K>::bounded || K == 1;
  const size_t d = test.dimension;
  const double* m0 = test.matrix.data();
  const double* m1 = m0 + d;
  double lo = 0.0, hi = 1.0;
  for (size_t i = 0; i < d && lo <= hi; ++i)
  {
    // (M d)_i = m0_i + w (m1_i - m0_i)
    const double k = m1[i] - m0[i];
    if (bounded)
    {
      if (k == 0.0) { if (m0[i] != 0.0) hi = -1.0; continue; }
      const double w = -m0[i] / k;
      lo = std::max(lo, w);
      hi = std::min(hi, w);
    }
    else if (k > 0.0) lo = std::max(lo, -m0[i] / k);
    else if (k < 0.0) hi = std::min(hi, -m0[i] / k);
    else if (m0[i] < 0.0) hi = -1.0;
  }
  if (!(lo <= hi)) return false;
  std::vector<double>& x = test.work_x;
  for (double w : {lo, hi})
  {
    const double slope = -(1.0 - w) * test.u[0] - w * test.u[1];
    if (!(slope > 0.0)) continue;
    for (double l = 1.0; l < 1e300; l *= 2.0)
    {
      x.assign({l * (1.0 - w), l * w});
      if (test.positive(x)) return true;
    }
  }
  return false;
}

template <class Model>
bool exact_two_constraints(const Decision<Model>& test)
{
  if (exact_one_constraint(test, 0) || exact_one_constraint(test, 1)) return true;

  const size_t d = test.dimension;
  const double* m0 = test.matrix.data();
  const double* m1 = m0 + d;
  std::vector<double>& x = test.work_x;
  x.assign(2, 0.0);
  if constexpr (dustib::model_id<Model> == 0)
  {
    // critical point: G x = h, G = M'M, h = -u - M'a
    double g00 = 0.0, g01 = 0.0, g11 = 0.0;
    double h0 = -test.u[0], h1 = -test.u[1];
    for (size_t i = 0; i < d; ++i)
    {
      g00 += m0[i] * m0[i];
      g01 += m0[i] * m1[i];
      g11 += m1[i] * m1[i];
      h0 -= m0[i] * test.a[i];
      h1 -= m1[i] * test.a[i];
    }
    const double det = g00 * g11 - g01 * g01;
    // dependent columns
    if (!(det > 1e-12 * g00 * g11)) return gauss_exact_search(test);
    x[0] = (g11 * h0 - g01 * h1) / det;
    x[1] = (g00 * h1 - g01 * h0) / det;
    return x[0] >= 0.0 && x[1] >= 0.0 && std::isfinite(x[0]) &&
      std::isfinite(x[1]) && test.positive(x);
  }
  else
  {
    if (d == 2)
    {
      // square system: M'y = -u, z = mean(y), M x = z - a
      const double det = m0[0] * m1[1] - m0[1] * m1[0];
      const double size = std::sqrt((m0[0] * m0[0] + m0[1] * m0[1]) *
                                    (m1[0] * m1[0] + m1[1] * m1[1]));
      if (std::isfinite(det) && std::abs(det) > 1e-12 * size)
      {
        const double y0 = (-test.u[0] * m1[1] + test.u[1] * m0[1]) / det;
        const double y1 = (-test.u[1] * m0[0] + test.u[0] * m1[0]) / det;
        const double z0 = Math<Model>::mean(y0), z1 = Math<Model>::mean(y1);
        if (Math<Model>::valid(z0) && Math<Model>::valid(z1))
        {
          const double r0 = z0 - test.a[0], r1 = z1 - test.a[1];
          x[0] = (r0 * m1[1] - m1[0] * r1) / det;
          x[1] = (m0[0] * r1 - r0 * m0[1]) / det;
          if (x[0] >= 0.0 && x[1] >= 0.0 && std::isfinite(x[0]) &&
              std::isfinite(x[1]) && test.positive(x)) return true;
        }
      }
    }
    return exact_two_unbounded(test);
  }
}

template <class Model>
bool exact_search(const Decision<Model>& test)
{
  if (test.constraints == 1) return exact_one_constraint(test, 0);
  if (test.constraints == 2) return exact_two_constraints(test);
  if constexpr (dustib::model_id<Model> == 0)
    return gauss_exact_search(test);
  else
    return false; // (at most 2 constraints here)
}

} // namespace dust_md

#endif
