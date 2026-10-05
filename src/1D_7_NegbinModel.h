#ifndef Negbin_1D_H
#define Negbin_1D_H

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

struct NegbinPolicy
{
  static inline double statistic(double data) { return data; }

  // See 1D_1_GaussModel.h for why this takes (point, a) rather than
  // (cumsum, point, t, s).
  static inline double costEval(double point, double a)
  {
    const double logTerm = -std::log1p(-std::exp(point));
    return logTerm - point * a;
  }

  static inline double costMin(const std::vector<double>& cumsum, unsigned int t, unsigned int s)
  {
    double res = 0;
    double m = (cumsum[t] - cumsum[s]) / (t - s);
    if(m > 0)
    {
      res = double(t - s) * std::log(1 + m) - (cumsum[t] - cumsum[s]) * std::log(m / (1 + m));
    }
    return res;
  }

  static inline double muMax(double a, double b)
  {
    double res = 1;
    if(b != 0){res = std::min(1.0, a / b);}
    return res;
  }

  static inline double xMax(double a, double b)
  {
    if (a < b) return -a / (a - b);
    return std::numeric_limits<double>::infinity();
  }

  // Tolerance-based: see 1D_5_BernModel.h for why exact equality isn't
  // reliable here (a reaches 0 via division, not just via summed zeros).
  static inline bool isLeftBoundary(double a) { return a < 1e-9; }
  static inline bool isRightBoundary(double a) { return false; } // unbounded above: no right boundary
  static inline double Dstar_leftboundary() { return 0; }
  static inline double Dstar_rightboundary() { return 0; } // unused: isRightBoundary always false
  static inline double Dstar_superLinearLimit() { return 0; }

  static inline double Dstar(double x)
  {
    const double inv1px = 1.0 / (1.0 + x);
    return x * std::log1p(-inv1px) - std::log1p(x);
  }
  static inline double DstarPrime(double x)
  {
    const double inv1px = 1.0 / (1.0 + x);
    return std::log1p(-inv1px);
  }
  static inline double DstarPrimeInv(double x) { return std::exp(x) / (1 - std::exp(x)); }
  static inline double DstarSecond(double x) { return 1.0 / x + 1.0 / (1.0 + x); }

  static inline const char* get_model() { return "negbin"; }
};

#endif
