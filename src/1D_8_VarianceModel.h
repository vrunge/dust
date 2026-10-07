#ifndef Variance_1D_H
#define Variance_1D_H

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

struct VariancePolicy
{
  /// statistic = data^2 (gaussian with mean 0)
  static inline double statistic(double data) { return data * data; }

  static inline double costEval(double point, double a)
  {
    return -0.5 * std::log(-2.0 * point) - point * a;
  }

  static inline double costMin(const std::vector<double>& cumsum, unsigned int t, unsigned int s)
  {
    double delta_t = t - s;
    double diff_cumsum = cumsum[t] - cumsum[s];
    if(diff_cumsum <= 0){diff_cumsum = 1e-100;} /// choice  1e-100 to avoid -Inf /// THIS IS IMPORTANT
    return 0.5 * delta_t * (1.0 + std::log(diff_cumsum / delta_t));
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

  // boundaries with a tolerance (a is obtained by a division)
  static inline bool isLeftBoundary(double a) { return a < 1e-9; }
  static inline bool isRightBoundary(double a) { return false; } // no right boundary
  static inline double Dstar_leftboundary() { return std::numeric_limits<double>::infinity(); }
  static inline double Dstar_rightboundary() { return 0; } // not used
  static inline double Dstar_superLinearLimit() { return 0; }

  static inline double Dstar(double x) { return -0.5 * (std::log(x) + 1.0); }
  static inline double DstarPrime(double x) { return -0.5 / x; }
  static inline double DstarPrimeInv(double x) { return -0.5 / x; }
  static inline double DstarSecond(double x) { return 0.5 / std::pow(x, 2); }

  static inline const char* get_model() { return "variance"; }
};

#endif
