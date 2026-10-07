/// Highway scan for MeanVar_DUST.cpp (compiled for each SIMD target)

#if defined(MEANVAR_DUST_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef MEANVAR_DUST_INL_H_
#undef MEANVAR_DUST_INL_H_
#else
#define MEANVAR_DUST_INL_H_
#endif

#include <hwy/highway.h>
#include <hwy/contrib/math/math-inl.h>

HWY_BEFORE_NAMESPACE();
namespace meanvar_hw {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

/// min of Q_s + cost(s, t) (from the last candidate: ties keep the largest s)
/// sums with rounding errors (*_lo)
std::pair<double, uint32_t> Scan(const double* pos, const double* sums, const double* sums_lo,
                                 const double* sums2, const double* sums2_lo,
                                 const double* costs, const uint32_t* active, size_t k,
                                 uint32_t t, double sum_t, double sum_t_lo,
                                 double sum2_t, double sum2_t_lo)
{
  const hn::ScalableTag<double> d;
  const size_t lanes = hn::Lanes(d);
  double best = inf;
  uint32_t index = 0;
  size_t i = k;
  double values[hn::MaxLanes(d)];
  while (i >= lanes) {
    i -= lanes;
    const auto len = hn::Sub(hn::Set(d, static_cast<double>(t)), hn::LoadU(d, pos + i));
    const auto s1 = hn::Add(hn::Sub(hn::Set(d, sum_t), hn::LoadU(d, sums + i)),
                            hn::Sub(hn::Set(d, sum_t_lo), hn::LoadU(d, sums_lo + i)));
    const auto s2 = hn::Add(hn::Sub(hn::Set(d, sum2_t), hn::LoadU(d, sums2 + i)),
                            hn::Sub(hn::Set(d, sum2_t_lo), hn::LoadU(d, sums2_lo + i)));
    const auto mean = hn::Div(s1, len);
    const auto variance = hn::Sub(hn::Div(s2, len), hn::Mul(mean, mean));
    const auto valid = hn::And(hn::Ge(len, hn::Set(d, 2.0)), hn::Gt(variance, hn::Zero(d)));
    const auto safe = hn::IfThenElse(valid, variance, hn::Set(d, 1.0));
    const auto candidate = hn::Add(hn::LoadU(d, costs + i),
      hn::Mul(hn::Set(d, 0.5), hn::Mul(len, hn::Add(hn::Set(d, 1.0), hn::CallLog(d, safe)))));
    hn::StoreU(hn::IfThenElse(valid, candidate, hn::Set(d, inf)), d, values);
    for (size_t j = lanes; j-- > 0;) {
      if (values[j] < best) { best = values[j]; index = active[i + j]; }
    }
  }
  while (i-- > 0) {
    const double value = costs[i] + segment_cost((sum_t - sums[i]) + (sum_t_lo - sums_lo[i]),
                                                 (sum2_t - sums2[i]) + (sum2_t_lo - sums2_lo[i]),
                                                 t - pos[i]);
    if (value < best) { best = value; index = active[i]; }
  }
  return {best, index};
}

} // namespace HWY_NAMESPACE
} // namespace meanvar_hw
HWY_AFTER_NAMESPACE();

#endif
