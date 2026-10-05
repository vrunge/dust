// Highway cost scan for MeanVar_DUST.cpp, compiled once per SIMD target
// (hwy/foreach_target.h); the best target is chosen at run time.

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

// Minimum over the k candidates of Q_s + segment_cost(s, t), scanning from
// the newest candidate so that ties keep the largest s.
std::pair<double, uint32_t> Scan(const double* pos, const double* sums, const double* sums2,
                                 const double* costs, const uint32_t* active, size_t k,
                                 uint32_t t, double sum_t, double sum2_t)
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
    const auto mean = hn::Div(hn::Sub(hn::Set(d, sum_t), hn::LoadU(d, sums + i)), len);
    const auto variance = hn::Sub(
      hn::Div(hn::Sub(hn::Set(d, sum2_t), hn::LoadU(d, sums2 + i)), len),
      hn::Mul(mean, mean));
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
    const double value = costs[i] + segment_cost(sum_t - sums[i], sum2_t - sums2[i], t - pos[i]);
    if (value < best) { best = value; index = active[i]; }
  }
  return {best, index};
}

} // namespace HWY_NAMESPACE
} // namespace meanvar_hw
HWY_AFTER_NAMESPACE();

#endif
