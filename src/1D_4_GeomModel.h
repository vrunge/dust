#ifndef Geom_1D_H
#define Geom_1D_H

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

struct GeomPolicy
{
  static inline double statistic(double data) { return data; }

  // See 1D_1_GaussModel.h for why this takes (point, a) rather than
  // (cumsum, point, t, s).
  static inline double costEval(double point, double a)
  {
    // Numerically stable computation of log(exp(-point) - 1), point < 0
    const double logTerm = -point + std::log1p(-std::exp(point));
    return -logTerm - point * a;
  }

  static inline double costMin(const std::vector<double>& cumsum, unsigned int t, unsigned int s)
  {
    double res = 0;
    double m = (cumsum[t] - cumsum[s]) / (t - s);
    if(m > 1)
    {
      res = (t - s) * std::log(m - 1) - (cumsum[t] - cumsum[s]) * std::log((m - 1) / m);
    }
    return res;
  }

  static inline double muMax(double a, double b)
  {
    if (b != 1) return std::min(1.0, (a - 1) / (b - 1));
    return 1.;
  }

  static inline double xMax(double a, double b)
  {
    if (a < b) return -(a - 1) / (a - b);
    return std::numeric_limits<double>::infinity();
  }

  // Tolerance-based: see 1D_5_BernModel.h for why exact equality isn't
  // reliable here (a reaches its boundary, 1, via division).
  static inline bool isLeftBoundary(double a) { return a < 1.0 + 1e-9; }
  static inline bool isRightBoundary(double a) { return false; } // unbounded above: no right boundary
  static inline double Dstar_leftboundary() { return 0; }
  static inline double Dstar_rightboundary() { return 0; } // unused: isRightBoundary always false

  static inline double Dstar(double x)
  {
    const double invx = 1.0 / x;
    return (x - 1.0) * std::log1p(-invx) - std::log(x);
  }
  static inline double DstarPrime(double x) { return std::log1p(-1.0 / x); }
  static inline double DstarPrimeInv(double x) { return 1 / (1 - std::exp(x)); }

  static inline const char* get_model() { return "geom"; }
};

#endif
