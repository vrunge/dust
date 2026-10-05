#ifndef Gauss_1D_H
#define Gauss_1D_H

#include <cmath>
#include <vector>
#include <string>
#include <limits>

// --------------------------------------------------------------------- //
// GaussPolicy: header-only, state-free cost-model policy for DUST_1D_T.
// Every member is `static inline`: no virtual dispatch, fully inlinable
// at every call site in the pruning loop (see 1D_DUST_Impl.h).
// Gaussian distribution with known (unit) variance.
// --------------------------------------------------------------------- //
struct GaussPolicy
{
  static inline double statistic(double data) { return data; }

  // a is the segment mean supplied by the pruning rule.
  static inline double costEval(double point, double a)
  {
    return point * (0.5 * point - a);
  }

  static inline double costMin(const std::vector<double>& cumsum, unsigned int t, unsigned int s)
  {
    return - 0.5 * (cumsum[t] - cumsum[s]) * (cumsum[t] - cumsum[s]) / (t - s);
  }

  static inline double muMax(double a, double b) { return 1; }
  static inline double xMax(double a, double b) { return std::numeric_limits<double>::infinity(); }

  static inline bool isLeftBoundary(double a) { return false; }
  static inline bool isRightBoundary(double a) { return false; } // unbounded domain: no right boundary
  static inline double Dstar_leftboundary() { return std::numeric_limits<double>::infinity(); }
  static inline double Dstar_rightboundary() { return std::numeric_limits<double>::infinity(); } // unused: isRightBoundary always false

  static inline double Dstar(double x) { return 0.5 * x * x; }
  static inline double DstarPrime(double x) { return x; }
  static inline double DstarPrimeInv(double x) { return x; }

  static inline const char* get_model() { return "gauss"; }
};

#endif
