#include <cmath>

#include "scan.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "dust/scan.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>
#include <hwy/contrib/math/math-inl.h>

HWY_BEFORE_NAMESPACE();
namespace dust {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

template <Cost K, class D, class V = hn::Vec<D>>
HWY_INLINE V SegmentCost(D d, V diff, V dt)
{
  const V zero = hn::Zero(d), one = hn::Set(d, 1.0);
  if constexpr (K == Cost::GaussianMean)
    return hn::Neg(hn::Div(hn::Mul(hn::Set(d, 0.5), hn::Mul(diff, diff)), dt));
  if constexpr (K == Cost::Poisson)
  {
    const V safe = hn::Max(diff, hn::Set(d, 1e-300));
    const V val = hn::Mul(diff, hn::Sub(one, hn::CallLog(d, hn::Div(safe, dt))));
    return hn::IfThenElse(hn::Gt(diff, zero), val, zero);
  }
  if constexpr (K == Cost::Exponential || K == Cost::GaussianVariance)
  {
    const V safe = hn::IfThenElse(hn::Gt(diff, zero), diff, hn::Set(d, 1e-100));
    const V val = hn::Mul(dt, hn::Add(one, hn::CallLog(d, hn::Div(safe, dt))));
    return K == Cost::GaussianVariance ? hn::Mul(hn::Set(d, 0.5), val) : val;
  }
  if constexpr (K == Cost::Geometric)
  {
    const V m = hn::Div(diff, dt);
    const V safe = hn::Max(hn::Sub(m, one), hn::Set(d, 1e-300));
    const V val = hn::Sub(hn::Mul(dt, hn::CallLog(d, safe)), hn::Mul(diff, hn::CallLog(d, hn::Div(safe, m))));
    return hn::IfThenElse(hn::Gt(m, one), val, zero);
  }
  if constexpr (K == Cost::Bernoulli || K == Cost::Binomial)
  {
    const V m = hn::Div(diff, dt);
    const V safe = hn::Max(hn::Min(m, hn::Set(d, 1.0 - 1e-300)), hn::Set(d, 1e-300));
    const V rest = hn::Sub(one, safe);
    const V val = hn::Neg(hn::Mul(dt, hn::Add(hn::Mul(safe, hn::CallLog(d, safe)), hn::Mul(rest, hn::CallLog(d, rest)))));
    return hn::IfThenElse(hn::And(hn::Gt(m, zero), hn::Lt(m, one)), val, zero);
  }
  if constexpr (K == Cost::NegativeBinomial)
  {
    const V m = hn::Div(diff, dt);
    const V safe = hn::Max(m, hn::Set(d, 1e-300));
    const V val = hn::Sub(hn::Mul(dt, hn::CallLog1p(d, safe)),
                          hn::Mul(diff, hn::CallLog(d, hn::Div(safe, hn::Add(one, safe)))));
    return hn::IfThenElse(hn::Gt(m, zero), val, zero);
  }
}

/// infinite for 1 point or a variance 0
template <class D, class V = hn::Vec<D>>
HWY_INLINE V MeanVarianceCost(D d, V sum, V sum2, V dt)
{
  const V mean = hn::Div(sum, dt);
  const V variance = hn::Sub(hn::Div(sum2, dt), hn::Mul(mean, mean));
  const auto valid = hn::And(hn::Ge(dt, hn::Set(d, 2.0)), hn::Gt(variance, hn::Zero(d)));
  const V log = hn::CallLog(d, hn::IfThenElse(valid, variance, hn::Set(d, 1.0)));
  const V val = hn::Mul(hn::Mul(hn::Set(d, 0.5), dt), hn::Add(hn::Set(d, 1.0), log));
  return hn::IfThenElse(valid, val, hn::Set(d, HUGE_VAL));
}

struct Arrays
{
  const double *pos, *cost, *hi, *lo;   // series row at hi + row * stride
  size_t stride, d;
};

/// m = lanes if kFull
template <Cost K, bool kConsecutive, bool kFull, class D, class V = hn::Vec<D>>
HWY_INLINE V CandidateCosts(D d, const Arrays& a, const double* hi_t, const double* lo_t,
                            V vt, V at, size_t i, size_t m)
{
  const auto load = [](D d, const double* p, size_t m) HWY_ATTR {
    if constexpr (kFull) return hn::LoadU(d, p);
    else return hn::LoadN(d, p, m);
  };
  // with the rounding errors of the sums of y^2
  const auto diff = [&](size_t row) HWY_ATTR {
    V sum = hn::Sub(hn::Set(d, hi_t[row]), load(d, a.hi + row * a.stride + i, m));
    if constexpr (K == Cost::GaussianVariance || K == Cost::GaussianMeanVariance)
      sum = hn::Add(sum, hn::Sub(hn::Set(d, lo_t[row]), load(d, a.lo + row * a.stride + i, m)));
    return sum;
  };
  const V dt = hn::Sub(vt, kConsecutive ? at : load(d, a.pos + i, m));
  V v = load(d, a.cost + i, m);
  for (size_t row = 0; row < a.d; row += Family<K>::statistics)
    if constexpr (K == Cost::GaussianMeanVariance) v = hn::Add(v, MeanVarianceCost(d, diff(row), diff(row + 1), dt));
    else v = hn::Add(v, SegmentCost<K>(d, diff(row), dt));
  return v;
}

template <Cost K, bool kConsecutive>
size_t ScanCandidates(const Candidates& c, const PrefixSums& sums, size_t t,
                      size_t b, size_t e, double* val)
{
  const hn::ScalableTag<double> d;
  const size_t N = hn::Lanes(d);
  const Arrays a{c.positions(), c.costs(), c.sums_hi(0), c.sums_lo(0), c.stride, c.d};
  const double *hi_t = sums.hi_at(t), *lo_t = sums.lo_at(t);
  const auto vt = hn::Set(d, static_cast<double>(t)), step = hn::Set(d, static_cast<double>(N));
  auto index = hn::Iota(d, static_cast<double>(b));
  auto at = hn::Iota(d, (c.size() ? a.pos[0] : 0.0) + static_cast<double>(b));
  auto best = hn::Set(d, HUGE_VAL);
  auto best_index = index;
  size_t i = b;
  for (; i + N <= e; i += N, index = hn::Add(index, step), at = hn::Add(at, step))
  {
    const auto v = CandidateCosts<K, kConsecutive, true>(d, a, hi_t, lo_t, vt, at, i, N);
    hn::StoreU(v, d, val + i);
    const auto better = hn::Lt(v, best);
    best = hn::IfThenElse(better, v, best);
    best_index = hn::IfThenElse(better, index, best_index);
  }
  if (i < e)
  {
    const auto v = CandidateCosts<K, kConsecutive, false>(d, a, hi_t, lo_t, vt, at, i, e - i);
    hn::StoreN(v, d, val + i, e - i);
    const auto better = hn::And(hn::Lt(v, best), hn::FirstN(d, e - i));
    best = hn::IfThenElse(better, v, best);
    best_index = hn::IfThenElse(better, index, best_index);
  }

  double vb[hn::MaxLanes(d)], ib[hn::MaxLanes(d)];
  hn::StoreU(best, d, vb);
  hn::StoreU(best_index, d, ib);
  size_t L = 0;
  for (size_t j = 1; j < N; j++)
    if (vb[j] < vb[L] || (vb[j] == vb[L] && ib[j] < ib[L])) L = j;
  return static_cast<size_t>(ib[L]);
}

ScanFunction GetScanFunction(Cost cost, bool consecutive)
{
  return with_cost(cost, [&](auto k) -> ScanFunction {
    constexpr Cost K = decltype(k)::value;
    return consecutive ? &ScanCandidates<K, true> : &ScanCandidates<K, false>;
  });
}

} // namespace HWY_NAMESPACE
} // namespace dust
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace dust {

HWY_EXPORT(GetScanFunction);

ScanFunction scan_function(Cost cost, bool consecutive)
{
  return HWY_DYNAMIC_DISPATCH(GetScanFunction)(cost, consecutive);
}

} // namespace dust
#endif
