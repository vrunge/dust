// Highway kernels for DUST_1D_HW.cpp. This file is compiled once per SIMD
// target (hwy/foreach_target.h) and the best one is chosen at run time.
//
// Only the O(k) optimal-cost scan is vectorized: its values also decide
// the PELT rule (Q_s + C(s,t) > Q_t) and the smallest-index rule. The
// DUST and DUSTib tests run scalar on the contiguous candidate arrays:
// the DUST active set stays small (about ten candidates), and evaluating
// every branch in every lane measured slower than the scalar early exits.

#if defined(DUST_1D_HW_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef DUST_1D_HW_INL_H_
#undef DUST_1D_HW_INL_H_
#else
#define DUST_1D_HW_INL_H_
#endif

#include <hwy/highway.h>
#include <hwy/contrib/math/math-inl.h>

HWY_BEFORE_NAMESPACE();
namespace dust_hw {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

// Optimal cost of a segment with statistic sum `diff` and length `dt`.
// K is dustib::model_id; the formulas match the scalar models' costMin.
template <int K, class D, class V = hn::Vec<D>>
HWY_INLINE V CostMin(D d, V diff, V dt)
{
  const V zero = hn::Zero(d), one = hn::Set(d, 1.0);
  if constexpr (K == 0) // gauss
    return hn::Neg(hn::Div(hn::Mul(hn::Set(d, 0.5), hn::Mul(diff, diff)), dt));
  if constexpr (K == 1) // poisson
  {
    const V safe = hn::Max(diff, hn::Set(d, 1e-300));
    const V val = hn::Mul(diff, hn::Sub(one, hn::CallLog(d, hn::Div(safe, dt))));
    return hn::IfThenElse(hn::Gt(diff, zero), val, zero);
  }
  if constexpr (K == 2 || K == 7) // exp, variance
  {
    const V safe = hn::IfThenElse(hn::Gt(diff, zero), diff, hn::Set(d, 1e-100));
    const V val = hn::Mul(dt, hn::Add(one, hn::CallLog(d, hn::Div(safe, dt))));
    return K == 7 ? hn::Mul(hn::Set(d, 0.5), val) : val;
  }
  if constexpr (K == 3) // geom
  {
    const V m = hn::Div(diff, dt);
    const V safe = hn::Max(hn::Sub(m, one), hn::Set(d, 1e-300));
    const V val = hn::Sub(hn::Mul(dt, hn::CallLog(d, safe)), hn::Mul(diff, hn::CallLog(d, hn::Div(safe, m))));
    return hn::IfThenElse(hn::Gt(m, one), val, zero);
  }
  if constexpr (K == 4 || K == 5) // bern, binom
  {
    const V m = hn::Div(diff, dt);
    const V safe = hn::Max(hn::Min(m, hn::Set(d, 1.0 - 1e-300)), hn::Set(d, 1e-300));
    const V rest = hn::Sub(one, safe);
    const V val = hn::Neg(hn::Mul(dt, hn::Add(hn::Mul(safe, hn::CallLog(d, safe)), hn::Mul(rest, hn::CallLog(d, rest)))));
    return hn::IfThenElse(hn::And(hn::Gt(m, zero), hn::Lt(m, one)), val, zero);
  }
  if constexpr (K == 6) // negbin
  {
    const V m = hn::Div(diff, dt);
    const V safe = hn::Max(m, hn::Set(d, 1e-300));
    const V val = hn::Sub(hn::Mul(dt, hn::CallLog1p(d, safe)), hn::Mul(diff, hn::CallLog(d, hn::Div(safe, hn::Add(one, safe)))));
    return hn::IfThenElse(hn::Gt(m, zero), val, zero);
  }
}

// Stores val[i] = Q_s + C(s, t) for every candidate s = A_pos[i] and
// returns the first index of the minimum.
template <int K>
size_t Scan(State& s, double ct, double tp)
{
  const hn::ScalableTag<double> d;
  const hn::RebindToSigned<decltype(d)> di;
  const size_t N = hn::Lanes(d), k = s.A_pos.size();
  const auto vct = hn::Set(d, ct), vtp = hn::Set(d, tp);
  auto best = hn::Set(d, HUGE_VAL);
  auto bestIdx = hn::Zero(di);
  for (size_t i = 0; i < k; i += N)
  {
    const size_t m = std::min(N, k - i);
    const auto v = hn::Add(hn::LoadN(d, s.A_cost.data() + i, m),
      CostMin<K>(d, hn::Sub(vct, hn::LoadN(d, s.A_cumsum.data() + i, m)),
                    hn::Sub(vtp, hn::LoadN(d, s.A_pos.data() + i, m))));
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

// Whether candidate i is pruned at time t with Q_t = qt.
template <int K, int M>
HWY_INLINE bool Prunes(const State& s, size_t i, double ct, double qt, double tp)
{
  // PELT rule: the oldest candidate has no constraint, and in PELT nothing has.
  if (M == PELT || i == 0) return s.val[i] > qt;
  const double n = tp - s.A_pos[i];
  const double a = (ct - s.A_cumsum[i]) / n, c = (qt - s.A_cost[i]) / n;
  if constexpr (M == DUST) return DualMax_DUST<ScalarModel<K>>::test_abcd(a, s.A_b[i], c, s.A_q[i]);
  else return dustib::test<K>(a, s.A_b[i], c, s.A_q[i]);
}

// Removes the dropped candidates and returns whether there were any.
// drop[i] then flags the survivors whose predecessor changed.
template <int M>
bool Compact(State& s)
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
      s.A_cumsum[w] = s.A_cumsum[i]; s.A_cost[w] = s.A_cost[i]; s.A_pos[w] = s.A_pos[i];
      if constexpr (M != PELT)
      {
        s.val[w] = s.val[i]; s.A_b[w] = s.A_b[i]; s.A_q[w] = s.A_q[i];
        if (prev && w > 0) s.set_constraint(w);
      }
      s.drop[w++] = prev;
    }
    prev = dropped;
  }
  s.A_cumsum.resize(w); s.A_cost.resize(w); s.A_pos.resize(w);
  if constexpr (M != PELT) { s.A_b.resize(w); s.A_q.resize(w); }
  return true;
}

template <int K, int M>
void Step(State& s, unsigned t)
{
  const double ct = s.cumsum[t], tp = t;
  const size_t k = s.A_pos.size();
  s.val.resize(k);
  const size_t arg = Scan<K>(s, ct, tp);
  const double qt = s.val[arg] + s.penalty;
  s.costRecord.push_back(qt);
  s.chptRecord.push_back(static_cast<int>(s.A_pos[arg]));

  if constexpr (M != OP)
  {
    s.drop.resize(k);
    for (size_t i = 0; i < k; i++) s.drop[i] = Prunes<K, M>(s, i, ct, qt, tp);
    // DUST iterates to a fixed point. A test depends only on the candidate
    // and its predecessor, so only the flagged survivors are re-tested.
    while (Compact<M>(s) && M == DUST)
    {
      bool any = false;
      for (size_t i = 0; i < s.A_pos.size(); i++)
        if (s.drop[i]) any |= (s.drop[i] = Prunes<K, M>(s, i, ct, qt, tp));
      if (!any) break;
    }
  }

  s.nb.push_back(static_cast<int>(s.A_pos.size()));
  s.push(ct, qt, tp, M == DUST || M == DUSTIB);
}

// Runs the dynamic program from the last processed time up to s.n.
void Advance(State& s)
{
  with_model(s.model, [&](auto model) {
    constexpr int K = decltype(model)::value;
    const unsigned from = static_cast<unsigned>(s.A_pos.back()) + 1;
    switch (s.method)
    {
      case DUST:   for (unsigned t = from; t <= s.n; t++) Step<K, DUST>(s, t);   break;
      case PELT:   for (unsigned t = from; t <= s.n; t++) Step<K, PELT>(s, t);   break;
      case OP:     for (unsigned t = from; t <= s.n; t++) Step<K, OP>(s, t);     break;
      case DUSTIB: for (unsigned t = from; t <= s.n; t++) Step<K, DUSTIB>(s, t); break;
    }
  });
}

} // namespace HWY_NAMESPACE
} // namespace dust_hw
HWY_AFTER_NAMESPACE();

#endif
