#ifndef DUST_COMMON_H
#define DUST_COMMON_H

#include <cmath>
#include <limits>
#include <vector>

namespace dust {

/// value > 0 beyond the rounding errors of terms of total size `scale`
inline bool positive(double value, double scale)
{
  constexpr double guard = 64 * std::numeric_limits<double>::epsilon();
  return std::isfinite(value) && std::isfinite(scale) && value > guard * (1 + scale);
}

/// prefix sums of d series at t * d + row, rounding errors in lo (TwoSum: tiny y^2 after large ones)
struct PrefixSums
{
  size_t d = 0;
  std::vector<double> hi, lo;

  void reset(size_t dimension)
  {
    d = dimension;
    hi.assign(d, 0.0);
    lo.assign(d, 0.0);
  }

  void add(const double* x)
  {
    const size_t at = hi.size() - d;
    for (size_t row = 0; row < d; ++row)
    {
      const double a = hi[at + row], sum = a + x[row], b = sum - a;
      hi.push_back(sum);
      lo.push_back(lo[at + row] + ((a - (sum - b)) + (x[row] - b)));
    }
  }

  const double* hi_at(size_t t) const { return hi.data() + t * d; }
  const double* lo_at(size_t t) const { return lo.data() + t * d; }

  double diff(size_t t, size_t s, size_t row) const
  {
    return (hi[t * d + row] - hi[s * d + row]) + (lo[t * d + row] - lo[s * d + row]);
  }
};

} // namespace dust

#endif
