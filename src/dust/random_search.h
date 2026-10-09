/// Multivariate pruning test: random evaluations of D (method "randomEval")

#ifndef DUST_RANDOM_SEARCH_H
#define DUST_RANDOM_SEARCH_H

#include <functional>

#include "decision.h"

namespace dust {

/// uniform random numbers on (0, 1)
using Uniform = std::function<double()>;

template <Cost K>
bool random_search(const Decision<K>& test, unsigned int draws, const Uniform& uniform)
{
  std::vector<double> x(test.constraints, 0.0);
  if (test.positive(x)) return true;
  for (unsigned int draw = 0; draw < draws; ++draw)
  {
    double total = 0.0;
    for (double& value : x)
    {
      value = -std::log(std::max(uniform(), 1e-15));
      total += value;
    }
    for (double& value : x) value /= total;

    double gradient = 0.0;
    double curvature = 0.0;
    bool interior = true;
    for (size_t i = 0; i < test.dimension; ++i)
    {
      if (Family<K>::boundary(test.a[i])) { interior = false; break; }
      double direction = 0.0;
      for (size_t j = 0; j < test.constraints; ++j)
        direction += test.matrix[j * test.dimension + i] * x[j];
      gradient -= Family<K>::theta(test.a[i]) * direction;
      curvature += Family<K>::curvature(test.a[i]) * direction * direction;
    }
    for (size_t j = 0; j < test.constraints; ++j) gradient -= test.u[j] * x[j];
    double radius = interior && gradient > 0.0 && curvature > 0.0 &&
      std::isfinite(curvature) ? gradient / curvature : 1.0;
    radius *= std::exp((2.0 * uniform() - 1.0) * std::log(4.0));
    if (!std::isfinite(radius) || radius <= 0.0) radius = 1.0;
    for (double& value : x) value *= radius;
    for (unsigned int attempt = 0; attempt < 40; ++attempt)
    {
      double scale = 0.0;
      const double score = test.value(x, scale);
      if (std::isfinite(score))
      {
        if (positive(score, scale)) return true;
        break;
      }
      for (double& value : x) value *= 0.5;
    }
  }
  return false;
}

} // namespace dust

#endif
