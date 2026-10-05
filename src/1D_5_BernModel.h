#ifndef Bern_1D_H
#define Bern_1D_H

#include <cmath>
#include <vector>
#include <string>

struct BernPolicy
{
  static inline double statistic(double data) { return data; }

  // See 1D_1_GaussModel.h for why this takes (point, a) rather than
  // (cumsum, point, t, s).
  static inline double costEval(double point, double a)
  {
    double val = 0;
    if (point <= 0.0) { val = std::log1p(std::exp(point)); }
    else { val = point + std::log1p(std::exp(-point)); }
    return val - point * a;
  }

  static inline double costMin(const std::vector<double>& cumsum, unsigned int t, unsigned int s)
  {
    double res = 0;
    double m = (cumsum[t] - cumsum[s]) / (t - s);
    if(m > 0 && m < 1)
    {
      res = - double(t - s) * (m * std::log(m) + (1 - m) * std::log(1 - m));
    }
    return res;
  }

  static inline double muMax(double a, double b)
  {
    double res = 1;
    if(b != 0 && b != 1){res = std::min(a / b, (1 - a) / (1 - b));}
    else{
      if(b == 0){res = 1 - a;} else {res = a;}
    }
    return res;
  }

  static inline double xMax(double a, double b)
  {
    if (a < b) return -a / (a - b);
    return (1 - a) / (a - b);
  }

  // Tolerance-based, not exact equality: `a` reaches the boundary via
  // division (e.g. (cumsum[t]-cumsum[s])/(t-s) landing on 1), which can
  // miss bit-exact equality by a couple of ULPs while still being well
  // within the region where Dstar/DstarPrime are numerically degenerate
  // (Dstar(1) is a literal 0*(-inf) = NaN; DstarPrime(x) diverges as
  // x -> 1). Bern/Binom are the only two models bounded on both sides,
  // hence the only two that need isRightBoundary at all.
  static inline bool isLeftBoundary(double a) { return a < 1e-9; }
  static inline bool isRightBoundary(double a) { return a > 1.0 - 1e-9; }
  static inline double Dstar_leftboundary() { return 0; }
  static inline double Dstar_rightboundary() { return 0; } // symmetric: Dstar(1) == Dstar(0)
  static inline double Dstar_superLinearLimit() { return 0; }

  static inline double Dstar(double x) { return x * std::log(x) + (1.0 - x) * std::log(1.0 - x); }
  static inline double DstarPrime(double x) { return std::log(x) - std::log(1.0 - x); }
  static inline double DstarPrimeInv(double x) { return std::exp(x) / (1 + std::exp(x)); }
  static inline double DstarSecond(double x) { return 1.0 / x + 1.0 / (1.0 - x); }

  static inline const char* get_model() { return "bern"; }
};

#endif
