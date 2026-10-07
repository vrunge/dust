/// MD DUST: random evaluations (method "randomEval")

#ifndef MD_RANDOM_EVAL_H
#define MD_RANDOM_EVAL_H

#include "MD_Decision.h"

namespace dust_md {

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

} // namespace dust_md

#endif
