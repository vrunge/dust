#ifndef DUST_SCAN_H
#define DUST_SCAN_H

#include <algorithm>
#include <vector>

#include "common.h"
#include "costs.h"

namespace dust {

/// candidates in increasing order: position, Q_s, prefix sums (one row of stride values per series);
/// live ones from `first`
struct Candidates
{
  size_t d = 0, first = 0, stride = 0;
  std::vector<double> pos, cost, hi, lo;

  void reset(size_t dimension)
  {
    d = dimension;
    first = stride = 0;
    pos.clear();
    cost.clear();
    hi.clear();
    lo.clear();
  }

  size_t size() const { return pos.size() - first; }
  const double* positions() const { return pos.data() + first; }
  const double* costs() const { return cost.data() + first; }
  const double* sums_hi(size_t row) const { return hi.data() + row * stride + first; }
  const double* sums_lo(size_t row) const { return lo.data() + row * stride + first; }
  size_t position(size_t i) const { return static_cast<size_t>(positions()[i]); }
  double cost_at(size_t i) const { return costs()[i]; }

  double sum(size_t row, size_t j, size_t i) const
  {
    const double *h = sums_hi(row), *l = sums_lo(row);
    return (h[i] - h[j]) + (l[i] - l[j]);
  }

  void push(size_t t, double q, const PrefixSums& sums)
  {
    const size_t n = pos.size();
    if (n == stride) grow(std::max<size_t>(2 * stride, 1));
    pos.push_back(static_cast<double>(t));
    cost.push_back(q);
    for (size_t row = 0; row < d; ++row)
    {
      hi[row * stride + n] = sums.hi_at(t)[row];
      lo[row * stride + n] = sums.lo_at(t)[row];
    }
  }

  /// pruning in place: move(i, j < i), then resize
  void move(size_t i, size_t j)
  {
    i += first;
    j += first;
    pos[j] = pos[i];
    cost[j] = cost[i];
    for (size_t row = 0; row < d; ++row)
    {
      hi[row * stride + j] = hi[row * stride + i];
      lo[row * stride + j] = lo[row * stride + i];
    }
  }

  void resize(size_t kept)
  {
    pos.resize(first + kept);
    cost.resize(first + kept);
  }

  /// memory moved when half is unused
  void drop_front(size_t k)
  {
    first += k;
    if (2 * first < pos.size()) return;
    const size_t n = size();
    for (size_t row = 0; row < d; ++row)
    {
      std::copy_n(hi.begin() + row * stride + first, n, hi.begin() + row * stride);
      std::copy_n(lo.begin() + row * stride + first, n, lo.begin() + row * stride);
    }
    pos.erase(pos.begin(), pos.begin() + first);
    cost.erase(cost.begin(), cost.begin() + first);
    first = 0;
  }

private:
  void grow(size_t new_stride)
  {
    std::vector<double> h(d * new_stride), l(d * new_stride);
    for (size_t row = 0; row < d; ++row)
    {
      std::copy_n(hi.begin() + row * stride, pos.size(), h.begin() + row * new_stride);
      std::copy_n(lo.begin() + row * stride, pos.size(), l.begin() + row * new_stride);
    }
    hi.swap(h);
    lo.swap(l);
    stride = new_stride;
  }
};

/// val[i] = Q_s + C(s, t) for the candidates b <= i < e, returns the first argmin
using ScanFunction = size_t (*)(const Candidates&, const PrefixSums&, size_t t,
                                size_t b, size_t e, double* val);

/// for the best SIMD target of this CPU (consecutive: positions s, s + 1, ... not loaded)
ScanFunction scan_function(Cost cost, bool consecutive);

} // namespace dust

#endif
