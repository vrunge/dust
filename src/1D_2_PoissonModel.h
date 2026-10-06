#ifndef Poisson_1D_H
#define Poisson_1D_H

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

struct PoissonPolicy
{
  static inline double statistic(double data) { return data; }

  // See 1D_1_GaussModel.h for why this takes (point, a) rather than
  // (cumsum, point, t, s).
  static inline double costEval(double point, double a)
  {
    return std::exp(point) - point * a;
  }

  static inline double costMin(const std::vector<double>& cumsum, unsigned int t, unsigned int s)
  {
    double res = 0.0;
    double diff = cumsum[t] - cumsum[s];
    if(diff > 0.0) /// should be != 0, we use > for robustness
    {
      res = diff * (1.0 - std::log(diff / (t - s)));
    }
    return res;
  }

  static inline double muMax(double a, double b)
  {
    if (b != 0) return std::min(1., a / b);
    return 1.;
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
  static inline double Dstar_rightboundary() { return std::numeric_limits<double>::infinity(); } // unused: isRightBoundary always false
  static inline double Dstar_superLinearLimit() { return std::numeric_limits<double>::infinity(); }

  static inline double Dstar(double x) { return (x * (std::log(x) - 1.0)); }
  static inline double DstarPrime(double x) { return std::log(x); }
  static inline double DstarPrimeInv(double x) { return std::exp(x); }
  static inline double DstarSecond(double x) { return std::pow(x, -1); }

  static inline const char* get_model() { return "poisson"; }
};

#endif
