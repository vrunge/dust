#include <Rcpp.h>
#include <cmath>

#include "1D_1_GaussModel.h"

using namespace Rcpp;

Gauss_1D::Gauss_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops)
  : DUST_1D(dualmax_algo, constr_index, nbLoops) {}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::costEval(double point, unsigned int t, unsigned int s) const
{
  //return 0.5*(t-s)*point*point - point*(cumsum[t] - cumsum[s]);
  const double dt  = static_cast<double>(t) - static_cast<double>(s);
  const double sum = cumsum[t] - cumsum[s];

  // Algebraic rewrite:
  // 0.5 * dt * point^2 - point * sum
  // = point * (0.5 * dt * point - sum)
  return point * (0.5 * point - sum/dt);
}

double Gauss_1D::costMin(unsigned int t, unsigned int s) const
{
  return - 0.5 * (cumsum[t] - cumsum[s]) * (cumsum[t] - cumsum[s]) / (t - s);
}

double Gauss_1D::statistic(double& data) const
{return(data);}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::dualEval(double point, double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  return (costRecord[s] - minCost_t) / (t - s) + point * (costRecord[s] - costRecord[r]) / (s - r)
    - 0.5 * std::pow((cumsum[t] - cumsum[s]) / (t - s) - point * ((cumsum[s] - cumsum[r]) / (s - r)), 2) / (1 - point);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::dualMax(double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  // Compute the optimal point on which to evaluate the duality function

  double a = (cumsum[t] - cumsum[s])/ (t - s); // Sbar_st
  double b = (cumsum[s] - cumsum[r])/ (s - r); // Sbar_rs
  double absAmB = std::abs(a - b);
  double sqrtB2p2C = std::sqrt(b*b + 2*(costRecord[s] - costRecord[r])/ (s - r));

  // Case 1: mu* = 0
  // deduce the following condition from the formula for mu*
  if (absAmB >= sqrtB2p2C)
    return (costRecord[s] - minCost_t) / (t - s)  - 0.5 * a*a;

  // Case 2: mu* > 0
    return (costRecord[s] - minCost_t) / (t - s) - 0.5*a*a + 0.5 * (absAmB - sqrtB2p2C)*(absAmB - sqrtB2p2C);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Gauss_1D::muMax(double a, double b) const
{
  return 1;
}

double Gauss_1D::xMax(double a, double b) const
{
  return std::numeric_limits<double>::infinity();
}

bool Gauss_1D::isLeftBoundary(double a) const {return false;}
double Gauss_1D::Dstar_leftboundary() const {return std::numeric_limits<double>::infinity();}

double Gauss_1D::Dstar_superLinearLimit() const {return std::numeric_limits<double>::infinity();}

double Gauss_1D::Dstar(double x) const
{
  return 0.5 * x * x;
}


double Gauss_1D::DstarPrime(double x) const
{
  return x;
}

double Gauss_1D::DstarPrimeInv(double x) const
{
  return x;
}

double Gauss_1D::DstarSecond(double x) const
{
  return 1.0;
}

std::string Gauss_1D::get_model() const { return "gauss"; }
