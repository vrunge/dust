/// Highway, VARIANCE model: Scan, Prunes, Compact and Step of
/// DUST_1D_HW-inl.h with the rounding errors of the sums (see 1D_Variance.h)
/// (included by DUST_1D_HW-inl.h)

#if defined(DUST_1D_HW_VARIANCE_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef DUST_1D_HW_VARIANCE_INL_H_
#undef DUST_1D_HW_VARIANCE_INL_H_
#else
#define DUST_1D_HW_VARIANCE_INL_H_
#endif

#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();
namespace dust_hw {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

constexpr int kVariance = 7;   // variance in Models

////////////////////////////////////////////////////////////////////////////////
inline size_t ScanVariance(VarianceState& s, double ct, double ct_lo, double tp)
{
  const hn::ScalableTag<double> d;
  const hn::RebindToSigned<decltype(d)> di;
  const size_t N = hn::Lanes(d), k = s.A_pos.size();
  const auto vct = hn::Set(d, ct), vct_lo = hn::Set(d, ct_lo), vtp = hn::Set(d, tp);
  auto best = hn::Set(d, HUGE_VAL);
  auto bestIdx = hn::Zero(di);
  for (size_t i = 0; i < k; i += N)
  {
    const size_t m = std::min(N, k - i);
    const auto statistic = hn::Add(hn::Sub(vct, hn::LoadN(d, s.A_cumsum.data() + i, m)),
                                   hn::Sub(vct_lo, hn::LoadN(d, s.A_cumsum_lo.data() + i, m)));
    const auto v = hn::Add(hn::LoadN(d, s.A_cost.data() + i, m),
      CostMin<kVariance>(d, statistic, hn::Sub(vtp, hn::LoadN(d, s.A_pos.data() + i, m))));
    hn::StoreN(v, d, s.val.data() + i, m);
    const auto better = hn::And(hn::Lt(v, best), hn::FirstN(d, m));
    best = hn::IfThenElse(better, v, best);
    bestIdx = hn::IfThenElse(hn::RebindMask(di, better), hn::Iota(di, static_cast<int64_t>(i)), bestIdx);
  }
  double vb[hn::MaxLanes(d)];
  int64_t ib[hn::MaxLanes(d)];
  hn::StoreU(best, d, vb);
  hn::StoreU(bestIdx, di, ib);
  size_t L = 0;
  for (size_t j = 1; j < N; j++)
    if (vb[j] < vb[L] || (vb[j] == vb[L] && ib[j] < ib[L])) L = j;
  return static_cast<size_t>(ib[L]);
}

////////////////////////////////////////////////////////////////////////////////
template <int M>
HWY_INLINE bool PrunesVariance(const VarianceState& s, size_t i, double ct, double ct_lo,
                               double qt, double tp)
{
  if (M == PELT || i == 0) return s.val[i] > qt;
  const double n = tp - s.A_pos[i];
  const double a = ((ct - s.A_cumsum[i]) + (ct_lo - s.A_cumsum_lo[i])) / n;
  const double c = (qt - s.A_cost[i]) / n;
  if constexpr (M == DUST) return DualMax_DUST<ScalarModel<kVariance>>::test_abcd(a, s.A_b[i], c, s.A_q[i]);
  else return dustib::test<kVariance>(a, s.A_b[i], c, s.A_q[i]);
}

////////////////////////////////////////////////////////////////////////////////
template <int M>
bool CompactVariance(VarianceState& s)
{
  const size_t k = s.A_pos.size();
  size_t w = std::find(s.drop.begin(), s.drop.begin() + k, 1) - s.drop.begin();
  if (w == k) return false;
  bool prev = false;
  for (size_t i = w; i < k; i++)
  {
    const bool dropped = s.drop[i];
    if (!dropped)
    {
      s.A_cumsum[w] = s.A_cumsum[i]; s.A_cumsum_lo[w] = s.A_cumsum_lo[i];
      s.A_cost[w] = s.A_cost[i]; s.A_pos[w] = s.A_pos[i];
      if constexpr (M != PELT)
      {
        s.val[w] = s.val[i]; s.A_b[w] = s.A_b[i]; s.A_q[w] = s.A_q[i];
        if (prev && w > 0) s.set_constraint_variance(w);
      }
      s.drop[w++] = prev;
    }
    prev = dropped;
  }
  s.A_cumsum.resize(w); s.A_cumsum_lo.resize(w); s.A_cost.resize(w); s.A_pos.resize(w);
  if constexpr (M != PELT) { s.A_b.resize(w); s.A_q.resize(w); }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
template <int M>
void StepVariance(VarianceState& s, unsigned t)
{
  const double ct = s.cumsum[t], ct_lo = s.cumsum_lo[t], tp = t;
  const size_t k = s.A_pos.size();
  s.val.resize(k);
  const size_t arg = ScanVariance(s, ct, ct_lo, tp);
  const double qt = s.val[arg] + s.penalty;
  s.costRecord.push_back(qt);
  s.chptRecord.push_back(static_cast<int>(s.A_pos[arg]));

  if constexpr (M != OP)
  {
    s.drop.resize(k);
    for (size_t i = 0; i < k; i++) s.drop[i] = PrunesVariance<M>(s, i, ct, ct_lo, qt, tp);
    while (CompactVariance<M>(s) && M == DUST)
    {
      bool any = false;
      for (size_t i = 0; i < s.A_pos.size(); i++)
        if (s.drop[i]) any |= (s.drop[i] = PrunesVariance<M>(s, i, ct, ct_lo, qt, tp));
      if (!any) break;
    }
  }

  s.nb.push_back(static_cast<int>(s.A_pos.size()));
  s.push_variance(ct, ct_lo, qt, tp, M == DUST || M == DUSTIB);
}

////////////////////////////////////////////////////////////////////////////////
void AdvanceVariance(VarianceState& s)
{
  const unsigned from = static_cast<unsigned>(s.A_pos.back()) + 1;
  switch (s.method)
  {
    case DUST:   for (unsigned t = from; t <= s.n; t++) StepVariance<DUST>(s, t);   break;
    case PELT:   for (unsigned t = from; t <= s.n; t++) StepVariance<PELT>(s, t);   break;
    case OP:     for (unsigned t = from; t <= s.n; t++) StepVariance<OP>(s, t);     break;
    case DUSTIB: for (unsigned t = from; t <= s.n; t++) StepVariance<DUSTIB>(s, t); break;
  }
}

} // namespace HWY_NAMESPACE
} // namespace dust_hw
HWY_AFTER_NAMESPACE();

#endif
