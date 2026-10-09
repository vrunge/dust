/// Univariate pruning tests with one constraint (Runge, Truong and Querne
/// 2025, arXiv:2507.02467). Candidate s is pruned at time t if the dual
///   D(mu) = -(c - mu d) - (1 - mu) A*((a - mu b) / (1 - mu))
/// is positive for some mu in [0, mu_max], with r < s the constraint and
///   a = mean on (s, t], b = mean on (r, s],
///   c = (Q_t - Q_s) / (t - s), d = (Q_s - Q_r) / (s - r).

#ifndef DUST_DECISION_1D_H
#define DUST_DECISION_1D_H

#include "common.h"
#include "costs.h"

namespace dust {

////////////////////////////////////////////////////////////////////////////////
/// Solver::Exact: closed-form maximum of the dual
template <Cost K>
bool closed_form_prunes(double a, double b, double c, double d)
{
  using M = Family<K>;
  constexpr double inf = std::numeric_limits<double>::infinity();
  // means within 1e-9 of a boundary are on it (a is obtained by a division)
  const bool left = K != Cost::GaussianMean && a < M::lower + 1e-9;
  const bool right = M::bounded && a > 1 - 1e-9;
  // A* on the boundary (infinite for gauss, exp and variance)
  constexpr double bound = K == Cost::GaussianMean || K == Cost::Exponential || K == Cost::GaussianVariance ? inf : 0.0;

  double mu_max = 1;
  if constexpr (M::bounded)
    mu_max = b == 0 ? 1 - a : b == 1 ? a : std::min(a / b, (1 - a) / (1 - b));
  else if constexpr (K != Cost::GaussianMean)
    if (b != M::lower) mu_max = std::min(1.0, (a - M::lower) / (b - M::lower));

  // one point segment or same means: D at mu = 0 and mu = mu_max
  if (left || right || std::abs(a - b) < 1e-14)
  {
    const double conj = left || right ? bound : M::conjugate(a);
    if (-c - conj > 0) return true;
    return std::abs(a - b) < 1e-14 && -(c - mu_max * d) - (1 - mu_max) * conj > 0;
  }

  // derivative at 0 negative: max at 0 (PELT test)
  if (-(a - b) * M::theta(a) - (c - d) < 0) return -M::conjugate(a) - c > 0;

  double x_max = inf;   // mu / (1 - mu) at mu_max
  if constexpr (K != Cost::GaussianMean)
  {
    if (a < b) x_max = -(a - M::lower) / (a - b);
    else if (M::bounded) x_max = (1 - a) / (a - b);
  }
  const double theta_star = -(c - d) / (a - b);
  if (x_max == inf)
  {
    // derivative at +infinity positive: D unbounded
    if (-M::theta(inf) - (c - d) > 0) return true;
  }
  else
  {
    // critical point outside the domain: max at x_max
    const double at_max = a + x_max * (a - b);
    const double critical = M::mean(theta_star);
    if (a > b ? critical >= at_max : critical <= at_max)
      return -bound - (c + x_max * (c - d)) > 0;
  }
  return M::partition(theta_star) - theta_star * a > c;
}

////////////////////////////////////////////////////////////////////////////////
/// Solver::ExactIB (DUSTib): same test with explicit domain checks and rounding guards
template <Cost K>
bool closed_form_ib_prunes(double a, double b, double c, double d)
{
  using M = Family<K>;
  if (!M::valid(a) || !M::valid(b) || !std::isfinite(c) || !std::isfinite(d)) return false;
  const double f = M::conjugate(a);
  if (positive(-f - c, std::abs(f) + std::abs(c))) return true;   // PELT
  if (a == b) return positive(d - c, std::abs(c) + std::abs(d));
  if (M::boundary(a)) return false;
  const double delta = a - b, e = c - d;
  // e == 0 needed: H grows like log, G goes to 0
  if (M::negative && delta > 0 && e <= 0) return true;
  const double u = delta * M::theta(a);
  if (!positive(-u - e, std::abs(u) + std::abs(c) + std::abs(d))) return false;
  const double r = -e / delta;
  if (!std::isfinite(r) || (M::negative && r >= 0)) return false;
  const double ar = M::partition(r), ra = r * a;
  return positive(ar - ra - c, std::abs(ar) + std::abs(ra) + std::abs(c));
}

////////////////////////////////////////////////////////////////////////////////
/// GaussianMeanVariance (one feature): m and m2 the means of y and y^2 and
/// q = (Q_t - Q_s) / (t - s) on (s, t], the same on (r, s] for the constraints r
struct Moments
{
  double m, m2, q;
};

namespace detail {

/// one constraint: maximum of the dual (x = 0 first if zero_first)
template <bool zero_first>
bool one_constraint(const Moments& a, const Moments& b)
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
  if constexpr (zero_first) {
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

/// dual with two constraints at (x1, x2)
inline double dual(const Moments& a, const Moments& b1, const Moments& b2, double x1, double x2)
{
  if (!(x1 >= 0.0 && x2 >= 0.0) || !std::isfinite(x1) || !std::isfinite(x2)) return -std::numeric_limits<double>::infinity();
  const double mean = a.m + x1 * (a.m - b1.m) + x2 * (a.m - b2.m);
  const double second = a.m2 + x1 * (a.m2 - b1.m2) + x2 * (a.m2 - b2.m2);
  const double variance = second - mean * mean;
  if (!(variance > 0.0) || !std::isfinite(variance)) return -std::numeric_limits<double>::infinity();
  return 0.5 * (1.0 + std::log(variance)) - a.q
    - x1 * (a.q - b1.q) - x2 * (a.q - b2.q);
}

} // namespace detail

/// one constraint
inline bool mean_variance_prunes(const Moments& a, const Moments& b)
{
  return detail::one_constraint<true>(a, b);
}

/// two constraints: the one-constraint maxima, then the critical point
inline bool mean_variance_prunes(const Moments& a, const Moments& b1, const Moments& b2)
{
  if (detail::one_constraint<false>(a, b1) || detail::one_constraint<false>(a, b2)) return true;
  if (!std::isfinite(a.q) || !std::isfinite(b1.q) || !std::isfinite(b2.q)) return false;

  const double a0 = a.m2, a1 = a.m2 - b1.m2, a2 = a.m2 - b2.m2;
  const double b0 = a.m, b_1 = a.m - b1.m, b_2 = a.m - b2.m;
  const double c1 = a.q - b1.q, c2 = a.q - b2.q;
  const double denomA = 2.0 * (b_1 * c2 - b_2 * c1);
  if (!std::isfinite(denomA) || std::abs(denomA) < 1e-12)
    return detail::dual(a, b1, b2, 0.0, 0.0) > 0.0;
  const double A = (b_1 * a2 - b_2 * a1) / denomA;
  double y;
  if (std::abs(b_1) >= std::abs(b_2) && std::abs(b_1) > 1e-12)
    y = (a1 - 2.0 * c1 * A) / (2.0 * b_1);
  else if (std::abs(b_2) > 1e-12)
    y = (a2 - 2.0 * c2 * A) / (2.0 * b_2);
  else return detail::dual(a, b1, b2, 0.0, 0.0) > 0.0;

  const double k1 = A + y * y - a0;
  const double k2 = y - b0;
  const double det = a1 * b_2 - a2 * b_1;
  if (!std::isfinite(det) || std::abs(det) < 1e-12)
    return detail::dual(a, b1, b2, 0.0, 0.0) > 0.0;
  const double x1 = (k1 * b_2 - k2 * a2) / det;
  const double x2 = (a1 * k2 - b_1 * k1) / det;
  return detail::dual(a, b1, b2, x1, x2) > 0.0;
}

} // namespace dust

#endif
