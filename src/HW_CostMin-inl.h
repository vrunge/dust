// Segment costs shared by the 1D and MD Highway candidate scans.
// Included inside each Highway target namespace.
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
    const V first = hn::Mul(dt, hn::CallLog1p(d, safe));
    const V second = hn::Mul(diff, hn::CallLog(d, hn::Div(safe, hn::Add(one, safe))));
    const V val = hn::Sub(first, second);
    return hn::IfThenElse(hn::Gt(m, zero), val, zero);
  }
}
