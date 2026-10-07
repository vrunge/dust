#ifndef DUST_1D_DUALMAX_POLICIES_H
#define DUST_1D_DUALMAX_POLICIES_H

#include <vector>
#include <cmath>
#include <limits>
#include "1D_DUSTib.h"

/// DUST pruning tests (Runge, Truong and Querne 2025, arXiv:2507.02467)

template <class Model>
inline bool isOnePointOrLinear(double a, double b)
{
  if(Model::isLeftBoundary(a) == true){return true;}
  if(Model::isRightBoundary(a) == true){return true;}
  /// the smallest visible gap is of size 1e-14
  if(std::abs(a - b) < 1e-14){return true;}
  return false;
}

// bern/binom bounded on both sides, the other models on the left (or not)
template <class Model>
inline bool specialCasePruning(double a,
                                double b,
                                double c,
                                double d,
                                double mu_max)
{
  if(Model::isLeftBoundary(a) == true)
  {
    if (-c  - Model::Dstar_leftboundary() > 0) {return true;} /// test in mu = 0
  }
  else if(Model::isRightBoundary(a) == true)
  {
    if (-c  - Model::Dstar_rightboundary() > 0) {return true;} /// test in mu = 0
  }
  if(std::abs(a - b) < 1e-14)
  {
    double Dstar_a = Model::isLeftBoundary(a) ? Model::Dstar_leftboundary()
                    : Model::isRightBoundary(a) ? Model::Dstar_rightboundary()
                    : Model::Dstar(a);
    if ( -c - Dstar_a > 0){return true;} /// mu = 0
    if (-(c - mu_max * d)  - (1 - mu_max)*Dstar_a > 0){return true;} /// mu = mu_max
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// DUAL = -(c - mu * d)  - (1- mu) * Dstar((a - mu * b) / (1 - mu));
//
// a = (cumsum[t] - cumsum[s]) / (t - s)
// b = (cumsum[s] - cumsum[r]) / (s - r)
// c = (Qt - Qs) / (t - s) =  (minCost_t - costRecord[s]) / (t - s);
// d = (Qs - Qr) / (s - r)  = (costRecord[s] - costRecord[r]) / (s - r);
////////////////////////////////////////////////////////////////////////////////

/// EXACT EVAL: closed-form maximum of the decision function.
template <class Model>
struct DualMax_DUST
{
  static constexpr const char* name() { return "DUST"; }

  static bool test_abcd(double a, double b, double c, double d)
  {
    double mu_max = Model::muMax(a, b);

    if(isOnePointOrLinear<Model>(a,b) == true){return specialCasePruning<Model>(a,b,c,d,mu_max);}

    // derivative in 0 negative: max in 0 (PELT)
    if(-(a-b)*Model::DstarPrime(a) - (c-d) < 0)
    {
      return (- Model::Dstar(a) -c > 0);
    }

    double x_max = Model::xMax(a,b);
    // DstarPrime(theta) = R
    double theta_star = -(c-d)/(a-b); // = R

    // derivative in +inf positive: dual unbounded, we prune
    if(x_max == std::numeric_limits<double>::infinity())
    {
      if(-Model::DstarPrime(std::numeric_limits<double>::infinity()) - (c-d) > 0)
      {
        return true;
      }
    }
    else
    {
      // critical point outside the domain?
      double theta_at_xmax = a + x_max*(a-b);
      double theta_crit = Model::DstarPrimeInv(theta_star);
      bool overshoot = (a > b) ? (theta_crit >= theta_at_xmax) : (theta_crit <= theta_at_xmax);
      if(overshoot)
      {
        double Dstar_at_xmax = (a > b) ? Model::Dstar_rightboundary() : Model::Dstar_leftboundary();
        return (-Dstar_at_xmax - (c + x_max*(c-d))) > 0;
      }
    }
    return(Model::costEval(theta_star, a) > c);
  }

  static bool test(const std::vector<double>& cumsum, const std::vector<double>& costRecord,
                    double minCost_t, unsigned int t, unsigned int s, unsigned int r)
  {
    double a = (cumsum[t] - cumsum[s]) / (t - s);
    double b = (cumsum[s] - cumsum[r]) / (s - r);
    double c = (minCost_t - costRecord[s]) / (t - s);
    double d = (costRecord[s] - costRecord[r]) / (s - r);
    return test_abcd(a, b, c, d);
  }
};

/// PELT = dual eval in mu = 0.
template <class Model>
struct DualMax_PELT
{
  static constexpr const char* name() { return "PELT"; }

  static bool test(const std::vector<double>& cumsum, const std::vector<double>& costRecord,
                    double minCost_t, unsigned int t, unsigned int s, unsigned int r)
  {
    double a = (cumsum[t] - cumsum[s]) / (t - s);
    double c = (minCost_t - costRecord[s]) / (t - s);
    if (- Model::Dstar(a) - c > 0) {return true;}
    return false;
  }
};

/// OP: see 1D_OP_Impl.h

#endif
