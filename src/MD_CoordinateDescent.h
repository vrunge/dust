/// MD DUST: coordinate descent (method "coordinateDescent")

#ifndef MD_COORDINATE_DESCENT_H
#define MD_COORDINATE_DESCENT_H

#include "MD_Decision.h"

namespace dust_md {

////////////////////////////////////////////////////////////////////////////////
/// GAUSS: D(x) = -|z|^2/2 - c - u.x, z = a + Mx (quadratic)
/// closed form update x_j <- max(0, x_j + g_j / |M_j|^2), g_j = -u_j - M_j.z
inline bool gauss_coordinate_search(const Decision<GaussPolicy>& test,
                                    unsigned int sweeps, double epsilon)
{
  const size_t d = test.dimension, k = test.constraints;
  double scale = 0.0;
  double previous_score = test.value(std::vector<double>(k, 0.0), scale);
  if (dustib::positive(previous_score, scale)) return true;
  if (!std::isfinite(previous_score)) return false;

  std::vector<double>& x = test.work_x;
  std::vector<double>& z = test.work_z;
  std::vector<double>& norm = test.work_norm;
  x.assign(k, 0.0);
  z.assign(test.a.begin(), test.a.end());
  norm.assign(k, 0.0);
  for (size_t j = 0; j < k; ++j)
    for (size_t i = 0; i < d; ++i)
      norm[j] += test.matrix[j * d + i] * test.matrix[j * d + i];
  double squares = 0.0;
  for (double v : z) squares += v * v;
  double linear = 0.0;     // u.x
  double linear_abs = 0.0; // sum |u_j x_j|

  for (unsigned int sweep = 0; sweep < sweeps; ++sweep)
  {
    bool changed = false;
    for (size_t j = 0; j < k; ++j)
    {
      const double* column = test.matrix.data() + j * d;
      double gradient = -test.u[j];
      for (size_t i = 0; i < d; ++i) gradient -= column[i] * z[i];
      if (norm[j] == 0.0)
      {
        // D linear in x_j
        if (!(gradient > 0.0)) continue;
        std::vector<double> ray(x);
        for (double step = 1.0; std::isfinite(step); step *= 2.0)
        {
          ray[j] = x[j] + step;
          if (test.positive(ray)) return true;
        }
        continue;
      }
      const double target = std::max(0.0, x[j] + gradient / norm[j]);
      const double delta = target - x[j];
      if (!(delta != 0.0) || !std::isfinite(target)) continue;
      changed = true;
      squares = 0.0;
      for (size_t i = 0; i < d; ++i)
      {
        z[i] += delta * column[i];
        squares += z[i] * z[i];
      }
      linear_abs += std::abs(test.u[j] * target) - std::abs(test.u[j] * x[j]);
      linear += test.u[j] * delta;
      x[j] = target;
      const double score = -0.5 * squares - test.c - linear;
      if (dustib::positive(score, std::abs(test.c) + 0.5 * squares + linear_abs) &&
          test.positive(x))
        return true;
    }
    if (!changed) break;
    if (no_pruning_certificate(test, z.data())) return false;
    if (epsilon >= 0.0)
    {
      const double score = -0.5 * squares - test.c - linear;
      if (!std::isfinite(score) || score - previous_score <= epsilon) break;
      previous_score = score;
    }
  }
  return false;
}

template <class Model>
bool coordinate_search(const Decision<Model>& test, unsigned int sweeps,
                       double epsilon)
{
  if constexpr (dustib::model_id<Model> == 0)
    return gauss_coordinate_search(test, sweeps, epsilon);
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

      const auto line = coordinate_argmax(test, base, j, lower, upper);
      const double at = line.first;
      if (line.second)
      {
        // D unbounded along x_j: we go further until D > 0
        for (double v = std::max(1.0, lower + 1.0); v < 1e300; v *= 2.0)
        {
          x[j] = v;
          if (test.coordinate_positive(base, j, v, linear_base, scale_base)) return true;
        }
        x[j] = old;
        continue;
      }
      x[j] = at;
      double new_scale = 0.0;
      const double new_score = test.coordinate_value(
        base, j, at, linear_base, scale_base, new_scale);
      if (dustib::positive(new_score, new_scale)) return true;
      if (new_score < old_score) x[j] = old;
      else changed = changed || x[j] != old;
    }
    if (!changed) break;
    // mean at the current point (for the certificate)
    std::vector<double>& z = test.work_z;
    z.assign(test.a.begin(), test.a.end());
    for (size_t j = 0; j < test.constraints; ++j)
      for (size_t i = 0; i < test.dimension; ++i)
        z[i] += test.matrix[j * test.dimension + i] * x[j];
    if (no_pruning_certificate(test, z.data())) return false;
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

} // namespace dust_md

#endif
