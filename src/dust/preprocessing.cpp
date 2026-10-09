#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "preprocessing.h"

namespace dust {

namespace {

bool all_finite(const std::vector<double>& y)
{
  return std::all_of(y.begin(), y.end(), [](double z) { return std::isfinite(z); });
}

// x is reordered
double median(std::vector<double>& x)
{
  const size_t n = x.size();
  std::nth_element(x.begin(), x.begin() + n / 2, x.end());
  double m = x[n / 2];
  if (n % 2 == 0)
  {
    std::nth_element(x.begin(), x.begin() + n / 2 - 1, x.end());
    m = (m + x[n / 2 - 1]) / 2.0;
  }
  return m;
}

} // namespace

double sd_diff(const std::vector<double>& y, const std::string& method)
{
  if (method == "HALL")
  {
    if (!all_finite(y) || y.size() < 5)
      throw std::invalid_argument("y is not a numeric vector or length < 5 (the HALL method cannot be used)");
    const double weight[4] = {0.1942, 0.2809, 0.3832, -0.8582};
    double sum_squares = 0.0;
    for (size_t j = 0; j + 3 < y.size(); ++j)
    {
      double column = 0.0;
      for (size_t i = 0; i < 4; ++i) column += weight[i] * y[j + i];
      sum_squares += column * column;
    }
    return std::sqrt(sum_squares / (y.size() - 3));
  }
  if (method != "MAD" && method != "SD") throw std::invalid_argument("method must be one of HALL, MAD, or SD");
  if (!all_finite(y) || y.size() < 2)
    throw std::invalid_argument("y is not a numeric vector or length < 2 (the " + method + " method cannot be used)");
  const size_t n = y.size() - 1;
  std::vector<double> diff(n);
  for (size_t i = 0; i < n; ++i) diff[i] = (y[i + 1] - y[i]) / std::sqrt(2.0);
  if (method == "MAD")
  {
    std::vector<double> sorted = diff;
    const double center = median(sorted);
    for (double& z : diff) z = std::abs(z - center);
    return median(diff) * 1.4826;
  }
  const double mean = std::accumulate(diff.begin(), diff.end(), 0.0) / n;
  double sum_squares = 0.0;
  for (double z : diff) sum_squares += std::pow(z - mean, 2);
  return std::sqrt(sum_squares / (n - 1));
}

std::vector<double> normalize_1d(std::vector<double> y, const std::string& type, double size)
{
  if (y.empty() || !all_finite(y)) throw std::invalid_argument("y must be a nonempty finite numeric vector");
  const auto any = [&](auto predicate) { return std::any_of(y.begin(), y.end(), predicate); };
  const auto divide = [&](double divisor) {
    for (double& z : y) z /= divisor;
    return y;
  };
  const double mean = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
  if (type == "variance")
  {
    for (double& z : y) z -= mean;
    if (any([](double z) { return z == 0.0; }))
      throw std::invalid_argument("variance normalization produced a zero residual");
    return y;
  }
  if (type == "gauss")
  {
    const double sd = sd_diff(y, "HALL");
    if (!(sd > 0.0) || !std::isfinite(sd))
      throw std::invalid_argument("Gaussian normalization requires positive finite noise scale");
    return divide(sd);
  }
  if (type == "poisson")
  {
    if (any([](double z) { return z < 0; })) throw std::range_error("negative data not compatible with poisson model");
    return mean == 0.0 ? y : divide(mean);
  }
  if (type == "exp")
  {
    if (any([](double z) { return z <= 0; })) throw std::range_error("Exponential data must be strictly positive");
    return divide(mean);
  }
  if (type == "binom")
  {
    if (!std::isfinite(size) || size <= 0 || std::floor(size) != size)
      throw std::invalid_argument("Binomial normalization requires a positive integer size (number of trials)");
    if (any([&](double z) { return z < 0 || z > size || std::floor(z) != z; }))
      throw std::invalid_argument("Binomial observations must be integer counts between 0 and size");
    return divide(size);
  }
  if (type == "negbin")
  {
    if (!std::isfinite(size) || size <= 0)
      throw std::invalid_argument("Negative Binomial normalization requires a positive known size");
    if (any([](double z) { return z < 0 || std::floor(z) != z; }))
      throw std::invalid_argument("Negative Binomial observations must be nonnegative integer counts");
    return divide(size);
  }
  if (type == "geom")
  {
    if (any([](double z) { return z < 1.0; })) throw std::invalid_argument("Geometric observations must be at least one");
    return y;
  }
  if (type == "bern")
  {
    if (any([](double z) { return z < 0.0 || z > 1.0; })) throw std::invalid_argument("Bernoulli observations must be in [0,1]");
    return y;
  }
  throw std::invalid_argument("Unsupported type specified.");
}

} // namespace dust
