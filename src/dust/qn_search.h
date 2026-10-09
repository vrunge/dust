#ifndef DUST_QN_SEARCH_H
#define DUST_QN_SEARCH_H

#include "decision.h"

namespace dust {

/// BFGS + Armijo line search, projected gradient if the BFGS step fails

template <Cost K>
bool qn_search(const Decision<K>& test, unsigned int iterations,
               double epsilon)
{
  const size_t p = test.constraints;
  std::vector<double> x(p, 0.0), gradient(p), next_gradient(p), trial(p);
  std::vector<double> direction(p), delta(p), y(p), hy(p);
  // mean on the boundary: fix the multipliers pushing it outside
  std::vector<bool> fixed(p, false);
  for (size_t i = 0; i < test.dimension; ++i)
    if (Family<K>::boundary(test.a[i]))
    {
      const bool upper = Family<K>::bounded && test.a[i] == 1.0;
      for (size_t j = 0; j < p; ++j)
      {
        const double b = test.matrix[j * test.dimension + i];
        if (upper ? b > 0.0 : b < 0.0) fixed[j] = true;
      }
    }
  double scale = 0.0;
  double score = test.value(x, scale);
  if (positive(score, scale)) return true;
  if (!std::isfinite(score)) return false;
  double trace = 0.0;
  if (!test.derivatives(x, gradient, &trace)) return false;
  for (size_t j = 0; j < p; ++j) if (fixed[j]) gradient[j] = 0.0;
  double max_gradient = 1.0;
  for (double g : gradient) max_gradient = std::max(max_gradient, std::abs(g));
  double initial_step = trace > 0.0 && std::isfinite(trace) ?
    1.0 / trace : 1.0 / max_gradient;
  if (!std::isfinite(initial_step) || initial_step <= 0.0)
    initial_step = 1.0 / max_gradient;

  std::vector<double> inverse(p * p, 0.0);
  const auto reset_inverse = [&]() {
    std::fill(inverse.begin(), inverse.end(), 0.0);
    for (size_t j = 0; j < p; ++j) inverse[j * p + j] = initial_step;
  };
  reset_inverse();

  for (unsigned int iteration = 0; iteration < iterations; ++iteration)
  {
    bool accepted = false;
    double next_score = score;
    // BFGS, then projected gradient
    for (unsigned int attempt = 0; attempt < 2u; ++attempt)
    {
      const bool bfgs_direction = attempt == 0;
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
        if (positive(next_score, scale)) return true;
        if (!std::isfinite(next_score) || !std::isfinite(slope) ||
            !std::isfinite(norm2) || !(slope > 0.0)) continue;
        const double required_gain = 0.01 * slope;
        if (!std::isfinite(required_gain) || next_score < score ||
            next_score < score + required_gain) continue;
        if (!test.derivatives(trial, next_gradient)) continue;
        for (size_t j = 0; j < p; ++j) if (fixed[j]) next_gradient[j] = 0.0;
        accepted = true;
        break;
      }
      if (accepted) break;
      reset_inverse();
    }
    if (!accepted) return false;
    if (epsilon >= 0.0 && next_score - score <= epsilon) return false;

    // BFGS update
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
      // no curvature (D almost linear): longer step
      else for (double& h : inverse) h *= 2.0;
    }
    x.swap(trial);
    gradient.swap(next_gradient);
    score = next_score;
    test.mean(x, test.work_z);
    if (no_pruning_certificate(test, test.work_z.data())) return false;
  }
  return false;
}

} // namespace dust

#endif
