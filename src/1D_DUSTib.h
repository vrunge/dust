#ifndef DUST_DUSTIB_H
#define DUST_DUSTIB_H

#include <cmath>
#include <limits>
#include <type_traits>
#include <stdexcept>
#include "1D_1_GaussModel.h"
#include "1D_2_PoissonModel.h"
#include "1D_3_ExpModel.h"
#include "1D_4_GeomModel.h"
#include "1D_5_BernModel.h"
#include "1D_6_BinomModel.h"
#include "1D_7_NegbinModel.h"
#include "1D_8_VarianceModel.h"

// One-constraint sign certificates for the eight cost models. They use
// the DUST dual decision function (Runge, Truong and Querné, 2025,
// doi:10.48550/arXiv.2507.02467) with explicit domain checks.

namespace dustib {
constexpr double guard = 64 * std::numeric_limits<double>::epsilon();
template<class M> constexpr int model_id =
  std::is_same<M, GaussPolicy>::value ? 0 :
  std::is_same<M, PoissonPolicy>::value ? 1 :
  std::is_same<M, ExpPolicy>::value ? 2 :
  std::is_same<M, GeomPolicy>::value ? 3 :
  std::is_same<M, BernPolicy>::value ? 4 :
  std::is_same<M, BinomPolicy>::value ? 5 :
  std::is_same<M, NegbinPolicy>::value ? 6 : 7;

template<int K> struct Math {
  static constexpr bool bounded = K == 4 || K == 5;
  static constexpr bool negative = K == 2 || K == 3 || K == 6 || K == 7;
  static bool valid(double a) {
    if (!std::isfinite(a)) return false;
    if constexpr (K == 0) return true;
    if constexpr (bounded) return a >= 0 && a <= 1;
    if constexpr (K == 3) return a >= 1;
    if constexpr (K == 2 || K == 7) return a > 0;
    return a >= 0;
  }
  static bool boundary(double a) {
    if constexpr (K == 0 || K == 2 || K == 7) return false;
    if constexpr (bounded) return a == 0 || a == 1;
    if constexpr (K == 3) return a == 1;
    return a == 0;
  }
  static double conjugate(double a) {
    if (boundary(a)) return 0;
    if constexpr (K == 0) return .5 * a * a;
    if constexpr (K == 1) return a * (std::log(a) - 1);
    if constexpr (K == 2) return -std::log(a) - 1;
    if constexpr (K == 7) return -.5 * (std::log(a) + 1);
    if constexpr (bounded) return a * std::log(a) + (1-a) * std::log1p(-a);
    if constexpr (K == 3) return (a-1) * std::log1p(-1/a) - std::log(a);
    if constexpr (K == 6) {
      // The log1p(-1/(1+a)) form loses small positive a to rounding.
      const double l = a < 1 ? std::log(a) - std::log1p(a) : -std::log1p(1/a);
      return a * l - std::log1p(a);
    }
  }
  static double theta(double a) {
    if constexpr (K == 0) return a;
    if constexpr (K == 1) return std::log(a);
    if constexpr (K == 2) return -1/a;
    if constexpr (K == 7) return -.5/a;
    if constexpr (bounded) return std::log(a) - std::log1p(-a);
    if constexpr (K == 3) return std::log1p(-1/a);
    if constexpr (K == 6) return a < 1 ? std::log(a)-std::log1p(a) : -std::log1p(1/a);
  }
  static double mean(double r) {
    if constexpr (K == 0) return r;
    if constexpr (K == 1) return std::exp(r);
    if constexpr (K == 2) return -1/r;
    if constexpr (K == 7) return -.5/r;
    if constexpr (bounded) { double z=std::exp(-std::abs(r)); return r >= 0 ? 1/(1+z) : z/(1+z); }
    if constexpr (K == 3) return -1/std::expm1(r);
    if constexpr (K == 6) return std::exp(r)/(-std::expm1(r));
  }
  static double partition(double r) {
    if constexpr (K == 0) return .5*r*r;
    if constexpr (K == 1) return std::exp(r);
    if constexpr (K == 2) return -std::log(-r);
    if constexpr (K == 7) return -.5*(std::log(-r)+std::log(2.0));
    if constexpr (bounded) return std::max(r,0.0)+std::log1p(std::exp(-std::abs(r)));
    if constexpr (K == 3 || K == 6) {
      double l = r < -std::log(2.0) ? std::log1p(-std::exp(r)) : std::log(-std::expm1(r));
      return (K == 3 ? r : 0) - l;
    }
  }
  static void validate(double y) {
    const double s = K == 7 ? y*y : y;
    if (!valid(s)) throw std::domain_error(
      "observation outside the model domain (Exponential and Variance require strictly positive sufficient statistics; Binomial requires proportions in [0,1])");
  }
};

inline bool positive(double value, double scale) {
  return std::isfinite(value) && std::isfinite(scale) && value > guard*(1+scale);
}

template<int K, int V = 2>
inline bool test(double a, double b, double c, double d) {
  using M = Math<K>;
  if (!M::valid(a) || !M::valid(b) || !std::isfinite(c) || !std::isfinite(d)) return false;
  const double f = M::conjugate(a);
  const double h0 = -f-c;
  const bool pelt = positive(h0, std::abs(f)+std::abs(c));
  if (pelt) return true;
  if (a == b) return positive(d-c, std::abs(c)+std::abs(d));
  if (M::boundary(a)) return pelt;
  const double delta = a-b, e = c-d;
  if constexpr (M::negative) {
    // e==0 is essential: H grows logarithmically, although G tends to 0.
    if (delta > 0 && e <= 0) return true;
  }
  const double u = delta*M::theta(a);
  if (!positive(-u-e, std::abs(u)+std::abs(c)+std::abs(d))) return pelt;
  const double r = -e/delta;
  if (!std::isfinite(r)) return pelt;
  if constexpr (M::negative) { if (r >= 0) return pelt; }
  const auto inequality = [&]() {
    const double ar=M::partition(r), ra=r*a;
    return pelt || positive(ar-ra-c,std::abs(ar)+std::abs(ra)+std::abs(c));
  };
  if constexpr (V == 2) return inequality();
  else {
    const double m = M::mean(r);
    // Under/overflow or a rounded endpoint: use the finite conjugate identity.
    if (!M::valid(m) || M::boundary(m) || !std::isfinite(m)) return inequality();
    const double x = (m-a)/delta;
    if (!(x > 0) || !std::isfinite(x)) return inequality();
    if constexpr (K != 0) {
      const double lo = K == 3 ? 1.0 : 0.0;
      if (delta < 0 && x >= (lo-a)/delta) return inequality();
      if constexpr (K == 4 || K == 5)
        if (delta > 0 && x >= (1-a)/delta) return inequality();
    }
    const double fm = M::conjugate(m);
    if constexpr (V == 1) {
      const double xe = x*e;
      return pelt || positive(-fm-c-xe,std::abs(fm)+std::abs(c)+std::abs(xe));
    } else {
      const double w = 1/(1+x), mu = x*w;
      const double wf = w*fm, md = mu*d;
      // mu rounding to 1 destroys the sign certificate at the excluded end.
      if (!(mu < 1) || w == 0) return inequality();
      return pelt || positive(-wf-c+md,std::abs(wf)+std::abs(c)+std::abs(md));
    }
  }
}
} // namespace dustib

template<class Model, int Version> struct DualMax_IB {
  static constexpr const char* name() { return "DUSTib"; }
  static bool test(const std::vector<double>& sums, const std::vector<double>& costs,
                   double qt, unsigned t, unsigned s, unsigned r) {
    return dustib::test<dustib::model_id<Model>,Version>(
      (sums[t]-sums[s])/(t-s), (sums[s]-sums[r])/(s-r),
      (qt-costs[s])/(t-s), (costs[s]-costs[r])/(s-r));
  }
};
template<class P> struct IsIB : std::false_type {};
template<class M, int V> struct IsIB<DualMax_IB<M,V>> : std::true_type {};
template<class Model> using DualMax_DUSTib = DualMax_IB<Model,2>;
#endif
