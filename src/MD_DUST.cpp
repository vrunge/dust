#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <utility>
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
  std::vector<std::vector<double>> sums; // one contiguous array per dimension

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
};

} // namespace md_hw

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "MD_DUST-inl.h"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include "MD_DUST-inl.h"
namespace md_hw {
HWY_EXPORT(Scan);
}

#include "DUST_1D_HW_Models.h"
#endif

namespace {

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

// See Runge, Truong and Querné (2025, arXiv:2507.02467) for the DUST test.
// Here D(x) = -sum A*(a + Mx) - c - u.x, for x >= 0.
template <class Model>
struct Decision
{
  size_t dimension;
  size_t constraints;
  bool highway;
  std::vector<double> a;
  std::vector<double> matrix; // one column per earlier index
  std::vector<double> u;
  double c;
  mutable std::vector<double> scratch;

  double value(const std::vector<double>& x, double& scale) const
  {
    // A pruning witness must lie in the multiplier orthant.
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

  // At a boundary mean, any nonzero comparison coefficient points outside
  // the model domain, so its multiplier is fixed at zero.
  bool derivatives(const std::vector<double>& x,
                   const std::vector<bool>& fixed,
                   std::vector<double>& gradient, double* trace = nullptr) const
  {
    for (size_t j = 0; j < constraints; ++j)
      gradient[j] = fixed[j] ? 0.0 : -u[j];
    if (trace) *trace = 0.0;
    for (size_t i = 0; i < dimension; ++i)
    {
      double z = a[i];
      bool varying = false;
      double norm2 = 0.0;
      for (size_t j = 0; j < constraints; ++j)
      {
        const double b = matrix[j * dimension + i];
        z += b * x[j];
        if (!fixed[j] && b != 0.0)
        {
          varying = true;
          if (trace) norm2 += b * b;
        }
      }
      if (!Math<Model>::valid(z)) return false;
      if (!varying) continue;
      if (Math<Model>::boundary(z)) return false;
      const double theta = Math<Model>::theta(z);
      if (!std::isfinite(theta)) return false;
      if (trace) *trace += norm2 * Model::DstarSecond(z);
      for (size_t j = 0; j < constraints; ++j)
        if (!fixed[j]) gradient[j] -= matrix[j * dimension + i] * theta;
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
      // Clipping is used for the search only; the test uses the true value.
      if constexpr (dustib::model_id<Model> == 3)
        z = std::max(z, 1.0 + 1e-12);
      else if constexpr (dustib::model_id<Model> == 4 || dustib::model_id<Model> == 5)
        z = std::min(1.0 - 1e-12, std::max(1e-12, z));
      else if constexpr (dustib::model_id<Model> != 0)
        z = std::max(z, 1e-12);
      result -= matrix[j * dimension + i] * Math<Model>::theta(z);
    }
    return result;
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

template <class Model>
bool coordinate_search(const Decision<Model>& test, unsigned int sweeps,
                       double epsilon)
{
  std::vector<double> x(test.constraints, 0.0);
  std::vector<double> base(test.dimension);
  double scale = 0.0;
  double previous_score = test.value(x, scale);
  if (dustib::positive(previous_score, scale)) return true;
  for (unsigned int sweep = 0; sweep < sweeps; ++sweep)
  {
    bool changed = false;
    for (size_t j = 0; j < test.constraints; ++j)
    {
      const double old = x[j];
      test.base_without(x, j, base);
      double linear_base = test.c;
      double scale_base = std::abs(test.c);
      for (size_t k = 0; k < test.constraints; ++k)
        if (k != j)
        {
          linear_base += test.u[k] * x[k];
          scale_base += std::abs(test.u[k] * x[k]);
        }
      double old_scale = 0.0;
      const double old_score = test.coordinate_value(
        base, j, old, linear_base, scale_base, old_scale);
      const auto limits = test.coordinate_bounds(base, j);
      const double lower = limits.first;
      const double upper = limits.second;
      if (!(upper > lower)) continue;
      x[j] = lower;
      const double at_lower = test.coordinate_derivative(base, j, lower);
      if (!(at_lower > 0.0))
      {
        double new_scale = 0.0;
        const double new_score = test.coordinate_value(
          base, j, lower, linear_base, scale_base, new_scale);
        if (new_score < old_score) x[j] = old;
        else if (old != lower) changed = true;
        continue;
      }

      double lo = lower;
      double hi = std::min(upper, std::max(lower + 1.0, 2.0 * old));
      if (!(hi > lower)) { x[j] = old; continue; }
      x[j] = hi;
      if (test.coordinate_positive(base, j, hi, linear_base, scale_base)) return true;
      double slope = test.coordinate_derivative(base, j, hi);
      for (unsigned int k = 0; slope > 0.0 && !std::isfinite(upper) && k < 40; ++k)
      {
        lo = hi;
        hi *= 2.0;
        if (!std::isfinite(hi)) break;
        x[j] = hi;
        if (test.coordinate_positive(base, j, hi, linear_base, scale_base)) return true;
        slope = test.coordinate_derivative(base, j, hi);
      }
      if (slope > 0.0 && std::isfinite(upper))
      {
        x[j] = upper;
        if (test.coordinate_positive(base, j, upper, linear_base, scale_base)) return true;
      }
      else if (slope > 0.0)
      {
        x[j] = old;
        continue;
      }
      else
      {
        for (unsigned int k = 0; k < 24; ++k)
        {
          const double mid = lo + 0.5 * (hi - lo);
          x[j] = mid;
          if (test.coordinate_positive(base, j, mid, linear_base, scale_base)) return true;
          if (test.coordinate_derivative(base, j, mid) > 0.0) lo = mid;
          else hi = mid;
        }
        x[j] = lo + 0.5 * (hi - lo);
      }
      double new_scale = 0.0;
      const double new_score = test.coordinate_value(
        base, j, x[j], linear_base, scale_base, new_scale);
      if (new_score < old_score) x[j] = old;
      else changed = changed || x[j] != old;
    }
    if (!changed) break;
    if (epsilon >= 0.0)
    {
      const double score = test.value(x, scale);
      if (dustib::positive(score, scale)) return true;
      if (!std::isfinite(score) || score - previous_score <= epsilon)
        break;
      previous_score = score;
    }
  }
  return false;
}

// Maximize the DUST decision over x >= 0. Iterative uses projected gradient;
// QN uses inverse BFGS with an Armijo line search. Both check every trial
// against the model domain and prune only at a positive feasible value.
template <class Model>
bool iterative_search(const Decision<Model>& test, unsigned int iterations,
                      bool quasi_newton, double epsilon)
{
  const size_t p = test.constraints;
  std::vector<double> x(p, 0.0), gradient(p), next_gradient(p), trial(p);
  std::vector<double> direction(p), delta(p), y(p), hy(p);
  std::vector<bool> fixed(p, false);
  double scale = 0.0;
  double score = test.value(x, scale);
  if (dustib::positive(score, scale)) return true;
  if (!std::isfinite(score)) return false;
  for (size_t i = 0; i < test.dimension; ++i)
    if (Math<Model>::boundary(test.a[i]))
      for (size_t j = 0; j < p; ++j)
        if (test.matrix[j * test.dimension + i] != 0.0) fixed[j] = true;

  double trace = 0.0;
  if (!test.derivatives(x, fixed, gradient, &trace)) return false;
  double max_gradient = 1.0;
  for (double g : gradient) max_gradient = std::max(max_gradient, std::abs(g));
  double initial_step = trace > 0.0 && std::isfinite(trace) ?
    1.0 / trace : 1.0 / max_gradient;
  if (!std::isfinite(initial_step) || initial_step <= 0.0)
    initial_step = 1.0 / max_gradient;

  std::vector<double> inverse(quasi_newton ? p * p : 0, 0.0);
  const auto reset_inverse = [&]() {
    std::fill(inverse.begin(), inverse.end(), 0.0);
    for (size_t j = 0; j < p; ++j) inverse[j * p + j] = initial_step;
  };
  if (quasi_newton) reset_inverse();

  for (unsigned int iteration = 0; iteration < iterations; ++iteration)
  {
    bool accepted = false;
    double next_score = score;
    // If BFGS fails, try projected gradient before ending the iteration.
    for (unsigned int attempt = 0; attempt < (quasi_newton ? 2u : 1u); ++attempt)
    {
      const bool bfgs_direction = quasi_newton && attempt == 0;
      for (size_t j = 0; j < p; ++j)
      {
        direction[j] = 0.0;
        if (fixed[j] || (x[j] == 0.0 && gradient[j] <= 0.0)) continue;
        if (bfgs_direction)
        {
          for (size_t k = 0; k < p; ++k)
            if (!fixed[k] && (x[k] > 0.0 || gradient[k] > 0.0))
              direction[j] += inverse[j * p + k] * gradient[k];
          if (x[j] == 0.0) direction[j] = std::max(0.0, direction[j]);
        }
        else direction[j] = gradient[j];
      }
      double alpha = bfgs_direction ? 1.0 : initial_step;
      for (unsigned int backtrack = 0; backtrack < 60; ++backtrack, alpha *= 0.5)
      {
        double slope = 0.0;
        double norm2 = 0.0;
        bool moved = false;
        bool finite = true;
        for (size_t j = 0; j < p; ++j)
        {
          const double proposed = x[j] + alpha * direction[j];
          if (!std::isfinite(proposed)) { finite = false; break; }
          trial[j] = std::max(0.0, proposed);
          delta[j] = trial[j] - x[j];
          slope += gradient[j] * delta[j];
          norm2 += delta[j] * delta[j];
          moved = moved || trial[j] != x[j];
        }
        if (!finite) continue;
        if (!moved) break;
        next_score = test.value(trial, scale);
        if (dustib::positive(next_score, scale)) return true;
        if (!std::isfinite(next_score) || !std::isfinite(slope) ||
            !std::isfinite(norm2) || !(slope > 0.0)) continue;
        const double required_gain = quasi_newton ? 0.01 * slope :
          slope - (0.5 / alpha) * norm2;
        if (!std::isfinite(required_gain) || next_score < score ||
            next_score < score + required_gain) continue;
        // Boundary points may witness pruning but cannot be optimizer states.
        if (!test.derivatives(trial, fixed, next_gradient)) continue;
        accepted = true;
        break;
      }
      if (accepted) break;
      if (quasi_newton) reset_inverse();
    }
    if (!accepted) return false;
    if (epsilon >= 0.0 && next_score - score <= epsilon) return false;

    if (quasi_newton)
    {
      double sy = 0.0, ss = 0.0, yy = 0.0;
      for (size_t j = 0; j < p; ++j)
      {
        y[j] = gradient[j] - next_gradient[j]; // gradient difference of F
        sy += delta[j] * y[j];
        ss += delta[j] * delta[j];
        yy += y[j] * y[j];
      }
      const double curvature_guard = 1e-10 * std::sqrt(ss) * std::sqrt(yy);
      if (std::isfinite(sy) && std::isfinite(curvature_guard) && sy > curvature_guard)
      {
        double yhy = 0.0;
        for (size_t j = 0; j < p; ++j)
        {
          hy[j] = 0.0;
          for (size_t k = 0; k < p; ++k) hy[j] += inverse[j * p + k] * y[k];
          yhy += y[j] * hy[j];
        }
        const double factor = (1.0 + yhy / sy) / sy;
        bool valid = std::isfinite(factor);
        for (size_t j = 0; j < p && valid; ++j)
          for (size_t k = 0; k < p; ++k)
          {
            double& h = inverse[j * p + k];
            h += factor * delta[j] * delta[k] -
              (hy[j] * delta[k] + delta[j] * hy[k]) / sy;
            if (!std::isfinite(h) || (j == k && !(h > 0.0))) valid = false;
          }
        if (!valid) reset_inverse();
      }
      else reset_inverse();
    }
    x.swap(trial);
    gradient.swap(next_gradient);
    score = next_score;
  }
  return false;
}

template <class Model>
bool random_search(const Decision<Model>& test, unsigned int draws)
{
  RNGScope scope;
  std::vector<double> x(test.constraints, 0.0);
  if (test.positive(x)) return true;
  for (unsigned int draw = 0; draw < draws; ++draw)
  {
    double total = 0.0;
    for (double& value : x)
    {
      value = -std::log(std::max(R::runif(0.0, 1.0), 1e-15));
      total += value;
    }
    for (double& value : x) value /= total;

    double gradient = 0.0;
    double curvature = 0.0;
    bool interior = true;
    for (size_t i = 0; i < test.dimension; ++i)
    {
      if (Math<Model>::boundary(test.a[i])) { interior = false; break; }
      double direction = 0.0;
      for (size_t j = 0; j < test.constraints; ++j)
        direction += test.matrix[j * test.dimension + i] * x[j];
      gradient -= Math<Model>::theta(test.a[i]) * direction;
      curvature += Model::DstarSecond(test.a[i]) * direction * direction;
    }
    for (size_t j = 0; j < test.constraints; ++j) gradient -= test.u[j] * x[j];
    double radius = interior && gradient > 0.0 && curvature > 0.0 &&
      std::isfinite(curvature) ? gradient / curvature : 1.0;
    radius *= std::exp((2.0 * R::runif(0.0, 1.0) - 1.0) * std::log(4.0));
    if (!std::isfinite(radius) || radius <= 0.0) radius = 1.0;
    for (double& value : x) value *= radius;
    for (unsigned int attempt = 0; attempt < 40; ++attempt)
    {
      double scale = 0.0;
      const double score = test.value(x, scale);
      if (std::isfinite(score))
      {
        if (dustib::positive(score, scale)) return true;
        break;
      }
      for (double& value : x) value *= 0.5;
    }
  }
  return false;
}

// In the Gaussian model the decision is quadratic:
// D(x) = D(0) + h'x - x'Gx/2, h = -M'S-u, G = M'M.
// On each face, solve G_II x_I = h_I and check the KKT conditions. The
// final pruning test evaluates the original decision at a feasible point.
enum class GaussMaximum { finite, unbounded, unresolved };

struct GaussMaximumResult
{
  GaussMaximum kind = GaussMaximum::unresolved;
  std::vector<double> point;
  double value = -std::numeric_limits<double>::infinity();
};

GaussMaximumResult gauss_joint_maximum(const Decision<GaussPolicy>& test)
{
  const size_t p = test.constraints;
  const size_t d = test.dimension;
  const double eps = std::numeric_limits<double>::epsilon();
  std::vector<double> h(p), gram(p * p);
  for (size_t j = 0; j < p; ++j)
  {
    h[j] = -test.u[j];
    for (size_t row = 0; row < d; ++row)
      h[j] -= test.a[row] * test.matrix[j * d + row];
    for (size_t k = 0; k <= j; ++k)
    {
      double dot = 0.0;
      for (size_t row = 0; row < d; ++row)
        dot += test.matrix[j * d + row] * test.matrix[k * d + row];
      gram[j * p + k] = gram[k * p + j] = dot;
    }
  }

  GaussMaximumResult result;
  std::vector<size_t> face;
  std::vector<double> unconstrained;
  std::vector<double> face_point, face_gradient;
  bool face_solved = false;
  std::function<bool(size_t)> visit = [&](size_t next) -> bool {
    face_solved = false;
    const size_t k = face.size();
    std::vector<double> q(d * k, 0.0), r(k * k, 0.0), x(p, 0.0);
    for (size_t col = 0; col < k; ++col)
    {
      const size_t j = face[col];
      double column_norm2 = 0.0;
      for (size_t row = 0; row < d; ++row)
      {
        q[col * d + row] = test.matrix[j * d + row];
        column_norm2 += q[col * d + row] * q[col * d + row];
      }
      // Reorthogonalize to keep the face solve stable.
      for (unsigned int pass = 0; pass < 2; ++pass)
        for (size_t previous = 0; previous < col; ++previous)
        {
          double projection = 0.0;
          for (size_t row = 0; row < d; ++row)
            projection += q[previous * d + row] * q[col * d + row];
          r[previous * k + col] += projection;
          for (size_t row = 0; row < d; ++row)
            q[col * d + row] -= projection * q[previous * d + row];
        }
      double remainder2 = 0.0;
      for (size_t row = 0; row < d; ++row)
        remainder2 += q[col * d + row] * q[col * d + row];
      const double diagonal = std::sqrt(remainder2);
      const double rank_guard = 64.0 * eps *
        static_cast<double>(std::max(d, p)) * std::sqrt(column_norm2);
      if (!std::isfinite(diagonal) || diagonal <= rank_guard)
      {
        // Check whether a dependent face has a nonnegative null direction
        // along which the decision increases.
        if (col + 1 == k)
        {
          std::vector<double> ray(p, 0.0);
          ray[j] = 1.0;
          for (size_t back = col; back-- > 0; )
          {
            double coefficient = r[back * k + col];
            for (size_t later = back + 1; later < col; ++later)
              coefficient += r[back * k + later] * ray[face[later]];
            ray[face[back]] = -coefficient / r[back * k + back];
          }
          bool nonnegative = true;
          for (size_t index : face)
          {
            if (ray[index] < -64.0 * eps) nonnegative = false;
            if (ray[index] < 0.0) ray[index] = 0.0;
          }
          double slope = 0.0, slope_scale = 0.0;
          for (size_t index : face)
          {
            slope += h[index] * ray[index];
            slope_scale += std::abs(h[index] * ray[index]);
          }
          if (nonnegative && dustib::positive(slope, slope_scale))
          {
            double scale = 0.0;
            const double at_zero = test.value(x, scale);
            double multiplier = (std::abs(at_zero) + scale + 1.0) / slope;
            if (!std::isfinite(multiplier) || multiplier <= 0.0)
              multiplier = 1.0;
            for (unsigned int attempt = 0; attempt < 64; ++attempt)
            {
              for (size_t index : face) x[index] = multiplier * ray[index];
              double witness_scale = 0.0;
              const double witness = test.value(x, witness_scale);
              if (dustib::positive(witness, witness_scale))
              {
                result = {GaussMaximum::unbounded, x, witness};
                return true;
              }
              if (!std::isfinite(multiplier * 2.0)) break;
              multiplier *= 2.0;
            }
          }
        }
        return false;
      }
      r[col * k + col] = diagonal;
      for (size_t row = 0; row < d; ++row)
        q[col * d + row] /= diagonal;
    }

    // Solve R'R x_I = h_I by two triangular solves.
    std::vector<double> intermediate(k, 0.0);
    for (size_t col = 0; col < k; ++col)
    {
      double rhs = h[face[col]];
      for (size_t previous = 0; previous < col; ++previous)
        rhs -= r[previous * k + col] * intermediate[previous];
      intermediate[col] = rhs / r[col * k + col];
    }
    bool feasible = true;
    for (size_t col = k; col-- > 0; )
    {
      double rhs = intermediate[col];
      for (size_t later = col + 1; later < k; ++later)
        rhs -= r[col * k + later] * x[face[later]];
      x[face[col]] = rhs / r[col * k + col];
      if (!std::isfinite(x[face[col]]) || x[face[col]] < 0.0)
        feasible = false;
    }
    if (p > 3)
    {
      if (k == p) unconstrained = x;
      face_point = x;
      face_gradient = h;
      for (size_t j = 0; j < p; ++j)
        for (size_t index : face)
          face_gradient[j] -= gram[j * p + index] * x[index];
      face_solved = true;
    }
    if (feasible)
    {
      for (size_t j = 0; j < p; ++j)
      {
        double gradient = p > 3 ? face_gradient[j] : h[j];
        double gradient_scale = std::abs(h[j]);
        for (size_t index : face)
        {
          const double term = gram[j * p + index] * x[index];
          if (p <= 3) gradient -= term;
          gradient_scale += std::abs(term);
        }
        const double tolerance = 128.0 * eps * (1.0 + gradient_scale);
        const bool active = std::find(face.begin(), face.end(), j) != face.end();
        if (!std::isfinite(gradient) ||
            (active ? std::abs(gradient) > tolerance :
                      gradient > tolerance))
        {
          feasible = false;
          break;
        }
      }
    }
    if (feasible)
    {
      double scale = 0.0;
      const double maximum = test.value(x, scale);
      if (std::isfinite(maximum))
      {
        result = {GaussMaximum::finite, x, maximum};
        return true;
      }
    }
    for (size_t j = next; j < p; ++j)
    {
      face.push_back(j);
      const bool found = visit(j + 1);
      face.pop_back();
      if (found) return true;
    }
    return false;
  };
  // Try the full face first; if it fails the KKT check, inspect other faces.
  if (p > 1)
  {
    face.resize(p);
    std::iota(face.begin(), face.end(), 0);
    if (visit(p)) return result;
    face.clear();
    // Follow active faces suggested by the unconstrained critical point.
    // Exhaustive enumeration handles cycles and unresolved faces.
    if (p > 3 && unconstrained.size() == p)
    {
      for (size_t j = 0; j < p; ++j)
        if (std::isfinite(unconstrained[j]) && unconstrained[j] > 0.0)
          face.push_back(j);
      std::vector<std::vector<size_t>> seen;
      for (size_t step = 0; step < 3 * p + 3; ++step)
      {
        if (std::find(seen.begin(), seen.end(), face) != seen.end()) break;
        seen.push_back(face);
        if (visit(p)) return result;
        if (!face_solved) break;
        size_t remove = p;
        double most_negative = 0.0;
        for (size_t j : face)
          if (face_point[j] < most_negative)
          {
            most_negative = face_point[j];
            remove = j;
          }
        if (remove < p)
        {
          face.erase(std::find(face.begin(), face.end(), remove));
          continue;
        }
        size_t add = p;
        double most_positive = 0.0;
        for (size_t j = 0; j < p; ++j)
          if (std::find(face.begin(), face.end(), j) == face.end() &&
              face_gradient[j] > most_positive)
          {
            most_positive = face_gradient[j];
            add = j;
          }
        if (add == p) break;
        face.insert(std::lower_bound(face.begin(), face.end(), add), add);
      }
      face.clear();
    }
  }
  visit(0);
  return result;
}

bool gauss_exact_search(const Decision<GaussPolicy>& test)
{
  const auto maximum = gauss_joint_maximum(test);
  return maximum.kind != GaussMaximum::unresolved &&
    test.positive(maximum.point);
}

class MDBase
{
public:
  MDBase(std::string model, std::string method, std::string backend,
         int constraints, int nb_iterations, double epsilon)
    : model_(std::move(model)), method_(std::move(method)),
      constraints_requested_(constraints), nb_iterations_(nb_iterations),
      epsilon_(epsilon)
  {
#ifdef HAVE_HIGHWAY
    highway_ = backend == "highway";
#else
    (void)backend;
#endif
  }
  virtual ~MDBase() = default;
  virtual void append_data(const NumericMatrix& data, Nullable<double> penalty) = 0;
  virtual void update_partition() = 0;
  virtual List get_partition() = 0;
  virtual List get_info() = 0;
  List dust(const NumericMatrix& data, Nullable<double> penalty)
  {
    append_data(data, penalty);
    update_partition();
    return get_partition();
  }

protected:
  std::string model_;
  std::string method_;
  int constraints_requested_;
  int nb_iterations_;
  double epsilon_;
  bool highway_ = false;
};

template <class Model>
class MDCommon : public MDBase
{
public:
  using MDBase::MDBase;

  void append_data(const NumericMatrix& data, Nullable<double> supplied) override
  {
    const size_t rows = data.nrow();
    const size_t cols = data.ncol();
    if (rows == 0) stop("data must have at least one row");
    if (dimension_ != 0 && rows != dimension_)
      stop("the number of rows cannot change after the first append");
    if (constraints_requested_ > static_cast<int>(rows))
      stop("constraints must be between 1 and the number of rows");
    if (cols == 0) return;
    for (double value : data) Math<Model>::validate(value);

    const double new_penalty = supplied.isNull() ?
      (dimension_ == 0 ? 2.0 * rows * std::log(static_cast<double>(cols)) : penalty_) :
      as<double>(supplied);
    if (!std::isfinite(new_penalty) || new_penalty < 0.0)
      stop("penalty must be a finite nonnegative number");
    if (dimension_ != 0 && new_penalty != penalty_)
      stop("penalty cannot change after the first nonempty append");

    if (dimension_ == 0)
    {
      dimension_ = rows;
      constraints_ = constraints_requested_ == 0 ? rows : constraints_requested_;
      penalty_ = new_penalty;
      sums_.assign(rows, 0.0);
      costs_.push_back(-penalty_);
      previous_.push_back(0);
      active_.push_back(0);
#ifdef HAVE_HIGHWAY
      if (highway_)
      {
        highway_state_.model = dustib::model_id<Model>;
        highway_state_.sums.resize(rows);
        highway_state_.push(0, costs_[0], sums_.data());
      }
#endif
    }
    sums_.reserve((length_ + cols + 1) * dimension_);
    for (size_t col = 0; col < cols; ++col)
      for (size_t row = 0; row < dimension_; ++row)
        sums_.push_back(sums_[(length_ + col) * dimension_ + row] +
                        Model::statistic(data(row, col)));
    length_ += cols;
  }

  List get_partition() override
  {
    if (length_ == 0) stop("append data before requesting a partition");
    if (done_ != length_) stop("update_partition before requesting a partition");
    std::vector<unsigned int> changepoints{static_cast<unsigned int>(length_)};
    for (unsigned int at = previous_[length_]; at != 0; at = previous_[at])
      changepoints.push_back(at);
    std::reverse(changepoints.begin(), changepoints.end());
    std::vector<unsigned int> remaining(active_.rbegin(), active_.rend());
    return List::create(
      _["changepoints"] = changepoints,
      _["lastIndexSet"] = remaining,
      _["backend"] = highway_ ? "highway" : "scalar",
      _["nb"] = counts_,
      _["costQ"] = NumericVector(costs_.begin() + 1, costs_.end())
    );
  }

  List get_info() override
  {
    RObject epsilon_info = R_NilValue;
    if (epsilon_ >= 0.0) epsilon_info = wrap(epsilon_);
    return List::create(
      _["backend"] = highway_ ? "highway" : "scalar",
      _["data_length"] = length_,
      _["dimension"] = dimension_,
      _["current_penalty"] = penalty_,
      _["model"] = model_,
      _["pruning_algo"] = method_,
      _["constraints"] = constraints_,
      _["nbIterations"] = nb_iterations_,
      _["epsilon"] = epsilon_info
    );
  }

protected:
  std::pair<double, unsigned int> scan(size_t t,
                                       std::vector<double>* candidates) const
  {
#ifdef HAVE_HIGHWAY
    if (highway_)
    {
      const auto result = HWY_DYNAMIC_DISPATCH(md_hw::Scan)(
        highway_state_, sums_.data() + t * dimension_,
        static_cast<double>(t), candidates);
      return {result.first, active_[result.second]};
    }
#endif
    if (candidates) candidates->resize(active_.size());
    double best = std::numeric_limits<double>::infinity();
    unsigned int argmin = 0;
    for (size_t i = 0; i < active_.size(); ++i)
    {
      const unsigned int s = active_[i];
      const double candidate = costs_[s] + cost(t, s);
      if (candidates) (*candidates)[i] = candidate;
      if (candidate < best) { best = candidate; argmin = s; }
    }
    return {best, argmin};
  }

  void retain_candidates(const std::vector<unsigned int>& survivors)
  {
#ifdef HAVE_HIGHWAY
    if (highway_) highway_state_.compact(survivors);
#else
    (void)survivors;
#endif
  }

  void push_candidate(size_t t)
  {
#ifdef HAVE_HIGHWAY
    if (highway_)
      highway_state_.push(t, costs_[t], sums_.data() + t * dimension_);
#else
    (void)t;
#endif
  }

  double cost(size_t t, size_t s) const
  {
    const double span = static_cast<double>(t - s);
    double total = 0.0;
    for (size_t row = 0; row < dimension_; ++row)
    {
      double mean = (sums_[t * dimension_ + row] - sums_[s * dimension_ + row]) / span;
      if constexpr (dustib::model_id<Model> == 4 || dustib::model_id<Model> == 5)
        mean = std::clamp(mean, 0.0, 1.0);
      total -= span * Math<Model>::conjugate(mean);
    }
    return total;
  }

  Decision<Model> decision(size_t t, size_t s,
                           const std::vector<unsigned int>& earlier) const
  {
    Decision<Model> test;
    test.dimension = dimension_;
    test.constraints = earlier.size();
    test.highway = highway_;
    test.a.resize(dimension_);
    test.scratch.resize(dimension_);
    test.matrix.resize(dimension_ * earlier.size());
    test.u.resize(earlier.size());
    const double span = static_cast<double>(t - s);
    test.c = (costs_[t] - costs_[s]) / span;
    for (size_t row = 0; row < dimension_; ++row)
      test.a[row] = (sums_[t * dimension_ + row] - sums_[s * dimension_ + row]) / span;
    for (size_t j = 0; j < earlier.size(); ++j)
    {
      const size_t r = earlier[j];
      const double earlier_span = static_cast<double>(s - r);
      const double prior_cost = (costs_[s] - costs_[r]) / earlier_span;
      test.u[j] = test.c - prior_cost;
      for (size_t row = 0; row < dimension_; ++row)
      {
        const double prior_mean = (sums_[s * dimension_ + row] -
                                   sums_[r * dimension_ + row]) / earlier_span;
        test.matrix[j * dimension_ + row] = test.a[row] - prior_mean;
      }
    }
    return test;
  }

  size_t dimension_ = 0;
  size_t constraints_ = 0;
  size_t length_ = 0;
  size_t done_ = 0;
  double penalty_ = 0.0;
  std::vector<double> sums_;
  std::vector<double> costs_;
  std::vector<unsigned int> previous_;
  std::vector<unsigned int> active_; // increasing order
  std::vector<int> counts_;
#ifdef HAVE_HIGHWAY
  md_hw::State highway_state_;
#endif
};

template <class Model>
class MDPruned : public MDCommon<Model>
{
public:
  using MDCommon<Model>::MDCommon;

  void update_partition() override
  {
    for (size_t t = this->done_ + 1; t <= this->length_; ++t)
    {
      std::vector<double> candidates;
      const auto best = this->scan(t, &candidates);
      this->costs_.push_back(best.first + this->penalty_);
      this->previous_.push_back(best.second);

      std::vector<unsigned int> survivors;
      survivors.reserve(this->active_.size() + 1);
      for (size_t position = 0; position < this->active_.size(); ++position)
      {
        const unsigned int s = this->active_[position];
        const double candidate = candidates[position];
        const double scale = std::abs(this->costs_[s]) +
          std::abs(candidate - this->costs_[s]) + std::abs(this->costs_[t]);
        bool prune = dustib::positive(candidate - this->costs_[t], scale);
        // Use as many earlier candidates as are currently available.
        if (!prune && this->method_ != "PELT" && !survivors.empty())
        {
          const size_t count = std::min(survivors.size(), this->constraints_);
          std::vector<unsigned int> earlier(
            survivors.end() - count, survivors.end());
          const auto test = this->decision(t, s, earlier);
          if (this->method_ == "coordinateDescent")
            prune = coordinate_search(test, this->nb_iterations_,
                                      this->epsilon_);
          else if (this->method_ == "iterative" || this->method_ == "QN")
            prune = iterative_search(test, this->nb_iterations_,
                                     this->method_ == "QN", this->epsilon_);
          else if constexpr (dustib::model_id<Model> == 0)
          {
            if (this->method_ == "exact")
              prune = gauss_exact_search(test);
            else
              prune = random_search(test, this->nb_iterations_);
          }
          else
            prune = random_search(test, this->nb_iterations_);
        }
        if (!prune)
        {
          survivors.push_back(s);
        }
      }
      this->counts_.push_back(static_cast<int>(survivors.size()));
      this->retain_candidates(survivors);
      survivors.push_back(static_cast<unsigned int>(t));
      this->active_.swap(survivors);
      this->push_candidate(t);
      this->done_ = t;
    }
  }
};

// Optimal partitioning keeps every candidate.
template <class Model>
class MDOP : public MDCommon<Model>
{
public:
  using MDCommon<Model>::MDCommon;

  void update_partition() override
  {
    for (size_t t = this->done_ + 1; t <= this->length_; ++t)
    {
      const auto best = this->scan(t, nullptr);
      this->costs_.push_back(best.first + this->penalty_);
      this->previous_.push_back(best.second);
      this->counts_.push_back(static_cast<int>(this->active_.size()));
      this->active_.push_back(static_cast<unsigned int>(t));
      this->push_candidate(t);
      this->done_ = t;
    }
  }
};

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
      method != "iterative" && method != "QN" &&
      method != "PELT" && method != "OP") stop("invalid MD method");
  if (backend != "highway" && backend != "scalar") stop("invalid MD backend");
  if (constraints < 0) stop("constraints must be positive or NULL");
  if (nb_iterations < 1) stop("nbIterations must be a positive integer");
  if (epsilon != -1.0 && (!std::isfinite(epsilon) || epsilon < 0.0))
    stop("epsilon must be finite, nonnegative, or NULL");
  if (epsilon >= 0.0 && method == "randomEval")
    stop("epsilon is not available for randomEval");
  const std::string effective = method == "exact" && model != "gauss" ? "PELT" : method;
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

} // namespace

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
