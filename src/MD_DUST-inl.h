// Highway candidate scan for all MD methods, compiled for each SIMD target.
#if defined(MD_DUST_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef MD_DUST_INL_H_
#undef MD_DUST_INL_H_
#else
#define MD_DUST_INL_H_
#endif

#include <hwy/highway.h>
#include <hwy/contrib/math/math-inl.h>

HWY_BEFORE_NAMESPACE();
namespace md_hw {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#include "HW_CostMin-inl.h"

template <int K>
std::pair<double, size_t> ScanModel(const State& s, const double* current,
                                     double t, std::vector<double>* values)
{
  const hn::ScalableTag<double> d;
  const size_t lanes = hn::Lanes(d), k = s.pos.size();
  const auto vt = hn::Set(d, t);
  double best = std::numeric_limits<double>::infinity();
  size_t argmin = 0;
  double block[hn::MaxLanes(d)];
  if (values) values->resize(k);
  for (size_t i = 0; i < k; i += lanes)
  {
    const size_t count = std::min(lanes, k - i);
    const auto span = hn::Sub(vt, hn::LoadN(d, s.pos.data() + i, count));
    auto cost = hn::Zero(d);
    for (size_t row = 0; row < s.sums.size(); ++row)
      cost = hn::Add(cost, CostMin<K>(
        d, hn::Sub(hn::Set(d, current[row]),
                   hn::LoadN(d, s.sums[row].data() + i, count)), span));
    cost = hn::Add(cost, hn::LoadN(d, s.cost.data() + i, count));
    hn::StoreN(cost, d, block, count);
    for (size_t j = 0; j < count; ++j)
    {
      const double candidate = block[j];
      if (values) (*values)[i + j] = candidate;
      if (candidate < best) { best = candidate; argmin = i + j; }
    }
  }
  return {best, argmin};
}

std::pair<double, size_t> Scan(const State& s, const double* current,
                                double t, std::vector<double>* values)
{
  switch (s.model)
  {
    case 0: return ScanModel<0>(s, current, t, values);
    case 1: return ScanModel<1>(s, current, t, values);
    case 2: return ScanModel<2>(s, current, t, values);
    case 3: return ScanModel<3>(s, current, t, values);
    case 4: return ScanModel<4>(s, current, t, values);
    case 5: return ScanModel<5>(s, current, t, values);
    case 6: return ScanModel<6>(s, current, t, values);
    default: return ScanModel<7>(s, current, t, values);
  }
}

} // namespace HWY_NAMESPACE
} // namespace md_hw
HWY_AFTER_NAMESPACE();

#endif
