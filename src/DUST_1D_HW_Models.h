// Vectorized cost models for the Highway backend.

#ifndef DUST_1D_HW_MODELS_H
#define DUST_1D_HW_MODELS_H

#include <hwy/highway.h>
#include <hwy/contrib/math/math-inl.h>

namespace hn = hwy::HWY_NAMESPACE;

// Gaussian model
struct VecGaussPolicy
{
  static inline double statistic(double x) { return x; }

  template <class D, class V> static inline V vDstar(D, V a) { return hn::Mul(hn::Set(D(), 0.5), hn::Mul(a, a)); }
  template <class D, class V> static inline V vDstarPrime(D, V a) { return a; }
  template <class D, class V> static inline V vMuMax(D d, V, V) { return hn::Set(d, 1.0); }
  template <class D, class V> static inline V vXMax(D d, V, V) { return hn::Set(d, std::numeric_limits<double>::infinity()); }
  template <class D, class V> static inline typename hn::Mask<D> vIsLeftBoundary(D d, V) { return hn::MaskFalse(d); }
  template <class D, class V> static inline typename hn::Mask<D> vIsRightBoundary(D d, V) { return hn::MaskFalse(d); }
  static constexpr double DstarLeftBoundary()  { return std::numeric_limits<double>::infinity(); }
  static constexpr double DstarRightBoundary() { return std::numeric_limits<double>::infinity(); }
  static constexpr double DstarPrimeAtInfinity() { return std::numeric_limits<double>::infinity(); }

  template <class D, class V>
  static inline V vCostEval(D, V point, V mean) { return hn::Mul(point, hn::Sub(hn::Mul(hn::Set(D(),0.5), point), mean)); }

  template <class D, class V>
  static inline V vCostMin(D d, V diff, V dt)
  {
    return hn::Neg(hn::Div(hn::Mul(hn::Set(d,0.5), hn::Mul(diff,diff)), dt));
  }
};

// Poisson model
struct VecPoissonPolicy
{
  static inline double statistic(double x) { return x; }

  template <class D, class V> static inline V vDstar(D d, V x) { return hn::Mul(x, hn::Sub(hn::CallLog(d, x), hn::Set(d,1.0))); }
  template <class D, class V> static inline V vDstarPrime(D d, V x) { return hn::CallLog(d, x); }
  template <class D, class V>
  static inline V vMuMax(D d, V a, V b)
  {
    auto ratio = hn::Div(a, b);
    auto m = hn::Min(hn::Set(d,1.0), ratio);
    return hn::IfThenElse(hn::Eq(b, hn::Zero(d)), hn::Set(d,1.0), m);
  }
  template <class D, class V>
  static inline V vXMax(D d, V a, V b)
  {
    auto finite = hn::Div(hn::Neg(a), hn::Sub(a,b));
    return hn::IfThenElse(hn::Lt(a,b), finite, hn::Set(d, std::numeric_limits<double>::infinity()));
  }
  template <class D, class V> static inline typename hn::Mask<D> vIsLeftBoundary(D d, V a) { return hn::Lt(a, hn::Set(d,1e-9)); }
  template <class D, class V> static inline typename hn::Mask<D> vIsRightBoundary(D d, V) { return hn::MaskFalse(d); }
  static constexpr double DstarLeftBoundary()  { return 0.0; }
  static constexpr double DstarRightBoundary() { return std::numeric_limits<double>::infinity(); }
  static constexpr double DstarPrimeAtInfinity() { return std::numeric_limits<double>::infinity(); }

  template <class D, class V>
  static inline V vCostEval(D d, V point, V mean) { return hn::Sub(hn::CallExp(d, point), hn::Mul(point, mean)); }

  template <class D, class V>
  static inline V vCostMin(D d, V diff, V dt)
  {
    auto safe = hn::Max(diff, hn::Set(d, 1e-300));
    auto val = hn::Mul(diff, hn::Sub(hn::Set(d,1.0), hn::CallLog(d, hn::Div(safe, dt))));
    return hn::IfThenElse(hn::Gt(diff, hn::Zero(d)), val, hn::Zero(d));
  }
};

// Exponential model
struct VecExpPolicy
{
  static inline double statistic(double x) { return x; }

  template <class D, class V> static inline V vDstar(D d, V x) { return hn::Neg(hn::Add(hn::CallLog(d, x), hn::Set(d,1.0))); }
  template <class D, class V> static inline V vDstarPrime(D d, V x) { return hn::Neg(hn::Div(hn::Set(d,1.0), x)); }
  template <class D, class V>
  static inline V vMuMax(D d, V a, V b)
  {
    auto m = hn::Min(hn::Set(d,1.0), hn::Div(a,b));
    return hn::IfThenElse(hn::Eq(b, hn::Zero(d)), hn::Set(d,1.0), m);
  }
  template <class D, class V>
  static inline V vXMax(D d, V a, V b)
  {
    auto finite = hn::Div(hn::Neg(a), hn::Sub(a,b));
    return hn::IfThenElse(hn::Lt(a,b), finite, hn::Set(d, std::numeric_limits<double>::infinity()));
  }
  template <class D, class V> static inline typename hn::Mask<D> vIsLeftBoundary(D d, V a) { return hn::Lt(a, hn::Set(d,1e-9)); }
  template <class D, class V> static inline typename hn::Mask<D> vIsRightBoundary(D d, V) { return hn::MaskFalse(d); }
  static constexpr double DstarLeftBoundary()  { return std::numeric_limits<double>::infinity(); }
  static constexpr double DstarRightBoundary() { return 0.0; }
  static constexpr double DstarPrimeAtInfinity() { return 0.0; }

  template <class D, class V>
  static inline V vCostEval(D d, V point, V mean) { return hn::Sub(hn::Neg(hn::CallLog(d, hn::Neg(point))), hn::Mul(point, mean)); }

  template <class D, class V>
  static inline V vCostMin(D d, V diff, V dt)
  {
    auto safe = hn::IfThenElse(hn::Gt(diff, hn::Zero(d)), diff, hn::Set(d, 1e-100));
    return hn::Mul(dt, hn::Add(hn::Set(d,1.0), hn::CallLog(d, hn::Div(safe, dt))));
  }
};

// Geometric model
struct VecGeomPolicy
{
  static inline double statistic(double x) { return x; }

  template <class D, class V>
  static inline V vDstar(D d, V x)
  {
    auto invx = hn::Div(hn::Set(d,1.0), x);
    return hn::Sub(hn::Mul(hn::Sub(x, hn::Set(d,1.0)), hn::CallLog1p(d, hn::Neg(invx))), hn::CallLog(d, x));
  }
  template <class D, class V> static inline V vDstarPrime(D d, V x) { return hn::CallLog1p(d, hn::Neg(hn::Div(hn::Set(d,1.0), x))); }
  template <class D, class V>
  static inline V vMuMax(D d, V a, V b)
  {
    auto m = hn::Min(hn::Set(d,1.0), hn::Div(hn::Sub(a,hn::Set(d,1.0)), hn::Sub(b,hn::Set(d,1.0))));
    return hn::IfThenElse(hn::Eq(b, hn::Set(d,1.0)), hn::Set(d,1.0), m);
  }
  template <class D, class V>
  static inline V vXMax(D d, V a, V b)
  {
    auto finite = hn::Div(hn::Neg(hn::Sub(a,hn::Set(d,1.0))), hn::Sub(a,b));
    return hn::IfThenElse(hn::Lt(a,b), finite, hn::Set(d, std::numeric_limits<double>::infinity()));
  }
  template <class D, class V> static inline typename hn::Mask<D> vIsLeftBoundary(D d, V a) { return hn::Lt(a, hn::Set(d, 1.0+1e-9)); }
  template <class D, class V> static inline typename hn::Mask<D> vIsRightBoundary(D d, V) { return hn::MaskFalse(d); }
  static constexpr double DstarLeftBoundary()  { return 0.0; }
  static constexpr double DstarRightBoundary() { return 0.0; }
  static constexpr double DstarPrimeAtInfinity() { return 0.0; }

  template <class D, class V>
  static inline V vCostEval(D d, V point, V mean)
  {
    auto logTerm = hn::Sub(hn::Neg(point), hn::Neg(hn::CallLog1p(d, hn::Neg(hn::CallExp(d, point)))));
    // logTerm = -point + log1p(-exp(point))  (matches scalar exactly)
    logTerm = hn::Add(hn::Neg(point), hn::CallLog1p(d, hn::Neg(hn::CallExp(d, point))));
    return hn::Sub(hn::Neg(logTerm), hn::Mul(point, mean));
  }

  template <class D, class V>
  static inline V vCostMin(D d, V diff, V dt)
  {
    auto m = hn::Div(diff, dt);
    auto mMinus1 = hn::Sub(m, hn::Set(d,1.0));
    auto safeMM1 = hn::Max(mMinus1, hn::Set(d,1e-300));
    auto val = hn::Sub(hn::Mul(dt, hn::CallLog(d, safeMM1)), hn::Mul(diff, hn::CallLog(d, hn::Div(safeMM1, m))));
    return hn::IfThenElse(hn::Gt(m, hn::Set(d,1.0)), val, hn::Zero(d));
  }
};

// Bernoulli and Binomial models
struct VecBernPolicy
{
  static inline double statistic(double x) { return x; }

  template <class D, class V> static inline V vDstar(D d, V x) { return hn::Add(hn::Mul(x, hn::CallLog(d,x)), hn::Mul(hn::Sub(hn::Set(d,1.0),x), hn::CallLog(d, hn::Sub(hn::Set(d,1.0),x)))); }
  template <class D, class V> static inline V vDstarPrime(D d, V x) { return hn::Sub(hn::CallLog(d,x), hn::CallLog(d, hn::Sub(hn::Set(d,1.0),x))); }
  template <class D, class V>
  static inline V vMuMax(D d, V a, V b)
  {
    auto general = hn::Min(hn::Div(a,b), hn::Div(hn::Sub(hn::Set(d,1.0),a), hn::Sub(hn::Set(d,1.0),b)));
    auto bIs0 = hn::Sub(hn::Set(d,1.0), a);
    auto bIs1 = a;
    auto maskB0 = hn::Eq(b, hn::Zero(d));
    auto maskB1 = hn::Eq(b, hn::Set(d,1.0));
    auto res = hn::IfThenElse(maskB0, bIs0, hn::IfThenElse(maskB1, bIs1, general));
    return res;
  }
  template <class D, class V>
  static inline V vXMax(D d, V a, V b)
  {
    auto lt = hn::Div(hn::Neg(a), hn::Sub(a,b));
    auto ge = hn::Div(hn::Sub(hn::Set(d,1.0),a), hn::Sub(a,b));
    return hn::IfThenElse(hn::Lt(a,b), lt, ge);
  }
  template <class D, class V> static inline typename hn::Mask<D> vIsLeftBoundary(D d, V a) { return hn::Lt(a, hn::Set(d,1e-9)); }
  template <class D, class V> static inline typename hn::Mask<D> vIsRightBoundary(D d, V a) { return hn::Gt(a, hn::Set(d, 1.0-1e-9)); }
  static constexpr double DstarLeftBoundary()  { return 0.0; }
  static constexpr double DstarRightBoundary() { return 0.0; }
  static constexpr double DstarPrimeAtInfinity() { return 0.0; } // unreachable: xMax never infinite for this model

  template <class D, class V>
  static inline V vCostEval(D d, V point, V mean)
  {
    auto neg = hn::CallLog1p(d, hn::CallExp(d, point));
    auto pos = hn::Add(point, hn::CallLog1p(d, hn::CallExp(d, hn::Neg(point))));
    auto val = hn::IfThenElse(hn::Le(point, hn::Zero(d)), neg, pos);
    return hn::Sub(val, hn::Mul(point, mean));
  }

  template <class D, class V>
  static inline V vCostMin(D d, V diff, V dt)
  {
    auto m = hn::Div(diff, dt);
    auto safeM = hn::Max(hn::Min(m, hn::Set(d, 1.0 - 1e-300)), hn::Set(d, 1e-300));
    auto val = hn::Neg(hn::Mul(dt, hn::Add(hn::Mul(safeM, hn::CallLog(d,safeM)),
                                            hn::Mul(hn::Sub(hn::Set(d,1.0),safeM), hn::CallLog(d, hn::Sub(hn::Set(d,1.0),safeM))))));
    auto inRange = hn::And(hn::Gt(m, hn::Zero(d)), hn::Lt(m, hn::Set(d,1.0)));
    return hn::IfThenElse(inRange, val, hn::Zero(d));
  }
};

using VecBinomPolicy = VecBernPolicy; // same formulas for normalized Binomial observations

// Negative Binomial model
struct VecNegbinPolicy
{
  static inline double statistic(double x) { return x; }

  template <class D, class V>
  static inline V vDstar(D d, V x)
  {
    auto inv1px = hn::Div(hn::Set(d,1.0), hn::Add(hn::Set(d,1.0), x));
    return hn::Sub(hn::Mul(x, hn::CallLog1p(d, hn::Neg(inv1px))), hn::CallLog1p(d, x));
  }
  template <class D, class V>
  static inline V vDstarPrime(D d, V x)
  {
    auto inv1px = hn::Div(hn::Set(d,1.0), hn::Add(hn::Set(d,1.0), x));
    return hn::CallLog1p(d, hn::Neg(inv1px));
  }
  template <class D, class V>
  static inline V vMuMax(D d, V a, V b)
  {
    auto m = hn::Min(hn::Set(d,1.0), hn::Div(a,b));
    return hn::IfThenElse(hn::Eq(b, hn::Zero(d)), hn::Set(d,1.0), m);
  }
  template <class D, class V>
  static inline V vXMax(D d, V a, V b)
  {
    auto finite = hn::Div(hn::Neg(a), hn::Sub(a,b));
    return hn::IfThenElse(hn::Lt(a,b), finite, hn::Set(d, std::numeric_limits<double>::infinity()));
  }
  template <class D, class V> static inline typename hn::Mask<D> vIsLeftBoundary(D d, V a) { return hn::Lt(a, hn::Set(d,1e-9)); }
  template <class D, class V> static inline typename hn::Mask<D> vIsRightBoundary(D d, V) { return hn::MaskFalse(d); }
  static constexpr double DstarLeftBoundary()  { return 0.0; }
  static constexpr double DstarRightBoundary() { return 0.0; }
  static constexpr double DstarPrimeAtInfinity() { return 0.0; }

  template <class D, class V>
  static inline V vCostEval(D d, V point, V mean)
  {
    auto logTerm = hn::Neg(hn::CallLog1p(d, hn::Neg(hn::CallExp(d, point))));
    return hn::Sub(logTerm, hn::Mul(point, mean));
  }

  template <class D, class V>
  static inline V vCostMin(D d, V diff, V dt)
  {
    auto m = hn::Div(diff, dt);
    auto safeM = hn::Max(m, hn::Set(d,1e-300));
    auto val = hn::Sub(hn::Mul(dt, hn::CallLog1p(d, safeM)), hn::Mul(diff, hn::CallLog(d, hn::Div(safeM, hn::Add(hn::Set(d,1.0),safeM)))));
    return hn::IfThenElse(hn::Gt(m, hn::Zero(d)), val, hn::Zero(d));
  }
};

// Variance model
struct VecVariancePolicy
{
  static inline double statistic(double x) { return x * x; }

  template <class D, class V> static inline V vDstar(D d, V x) { return hn::Mul(hn::Set(d,-0.5), hn::Add(hn::CallLog(d,x), hn::Set(d,1.0))); }
  template <class D, class V> static inline V vDstarPrime(D d, V x) { return hn::Neg(hn::Div(hn::Set(d,0.5), x)); }
  template <class D, class V>
  static inline V vMuMax(D d, V a, V b)
  {
    auto m = hn::Min(hn::Set(d,1.0), hn::Div(a,b));
    return hn::IfThenElse(hn::Eq(b, hn::Zero(d)), hn::Set(d,1.0), m);
  }
  template <class D, class V>
  static inline V vXMax(D d, V a, V b)
  {
    auto finite = hn::Div(hn::Neg(a), hn::Sub(a,b));
    return hn::IfThenElse(hn::Lt(a,b), finite, hn::Set(d, std::numeric_limits<double>::infinity()));
  }
  template <class D, class V> static inline typename hn::Mask<D> vIsLeftBoundary(D d, V a) { return hn::Lt(a, hn::Set(d,1e-9)); }
  template <class D, class V> static inline typename hn::Mask<D> vIsRightBoundary(D d, V) { return hn::MaskFalse(d); }
  static constexpr double DstarLeftBoundary()  { return std::numeric_limits<double>::infinity(); }
  static constexpr double DstarRightBoundary() { return 0.0; }
  static constexpr double DstarPrimeAtInfinity() { return 0.0; }

  template <class D, class V>
  static inline V vCostEval(D d, V point, V mean) { return hn::Sub(hn::Mul(hn::Set(d,-0.5), hn::CallLog(d, hn::Mul(hn::Set(d,-2.0),point))), hn::Mul(point, mean)); }

  template <class D, class V>
  static inline V vCostMin(D d, V diff, V dt)
  {
    auto safe = hn::IfThenElse(hn::Gt(diff, hn::Zero(d)), diff, hn::Set(d, 1e-100));
    return hn::Mul(hn::Set(d,0.5), hn::Mul(dt, hn::Add(hn::Set(d,1.0), hn::CallLog(d, hn::Div(safe, dt)))));
  }
};

#endif
