#ifndef DUST_1D_DUALMAX_POLICIES_H
#define DUST_1D_DUALMAX_POLICIES_H

#include <vector>
#include <cmath>
#include <limits>
#include "1D_DUSTib.h"

// Pruning rules for the scalar engine. The DUST rule evaluates the dual
// decision function of Runge, Truong and Querné (2025),
// doi:10.48550/arXiv.2507.02467.

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

  // Shared by the scalar engine (test() below) and the Highway engine,
  // which computes a, b, c, d from its contiguous candidate arrays.
  static bool test_abcd(double a, double b, double c, double d)
  {
    // Boundary a (Dstar takes its limit value) or linear dual (a == b up to
    // the smallest visible gap 1e-14): test mu = 0, and mu = mu_max if linear.
    const bool left = Model::isLeftBoundary(a), right = !left && Model::isRightBoundary(a);
    const bool linear = std::abs(a - b) < 1e-14;
    if (left || right || linear)
    {
      const double Dstar_a = left ? Model::Dstar_leftboundary() : right ? Model::Dstar_rightboundary() : Model::Dstar(a);
      if (-c - Dstar_a > 0) return true;
      if (!linear) return false;
      const double mu_max = Model::muMax(a, b);
      return -(c - mu_max * d) - (1 - mu_max) * Dstar_a > 0;
    }

    // if derivative in 0 is negative, the constrained maximum is at x=0
    // (the PELT point): a is guaranteed away from either boundary here,
    // so DstarPrime(a)/Dstar(a) are both well defined.
    if(-(a-b)*Model::DstarPrime(a) - (c-d) < 0)
    {
      return (- Model::Dstar(a) -c > 0);
    }

    double x_max = Model::xMax(a,b);
    // theta_star stores the right-hand side of DstarPrime(theta) = R.
    // costEval takes R directly; the domain check needs DstarPrimeInv(R).
    double theta_star = -(c-d)/(a-b); // = R

    // if derivative in +inf is positive, the dual is unbounded above on
    // an infinite domain: pruning succeeds without needing theta_star.
    if(x_max == std::numeric_limits<double>::infinity())
    {
      if(-Model::DstarPrime(std::numeric_limits<double>::infinity()) - (c-d) > 0)
      {
        return true;
      }
    }
    else
    {
      // Check whether the stationary point lies beyond the model domain.
      double theta_at_xmax = a + x_max*(a-b);
      double theta_crit = Model::DstarPrimeInv(theta_star); // actual critical theta
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
    return -Model::Dstar(a) - c > 0;
  }
};

// OP uses DUST_1D_OP_T in 1D_OP_Impl.h so that no candidate is pruned.

#endif
