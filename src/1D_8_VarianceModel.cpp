#include <Rcpp.h>
#include <cmath>

#include "1D_8_VarianceModel.h"

using namespace Rcpp;

Variance_1D::Variance_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops)
  : DUST_1D(dualmax_algo, constr_index, nbLoops) {}


double Variance_1D::statistic(double& data) const
{return(data * data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Variance_1D::costEval(double point, unsigned int t, unsigned int s) const
{
  return -0.5*std::log(-2.0*point) - point*(cumsum[t] - cumsum[s])/(1.0*(t-s));
}

double Variance_1D::costMin(unsigned int t, unsigned int s) const
{
  double delta_t = t - s;
  double diff_cumsum = cumsum[t] - cumsum[s];
  if(diff_cumsum <= 0){diff_cumsum = 1e-100;} /// choice  1e-100 to avoid -Inf /// THIS IS IMPORTANT
  return 0.5 * delta_t * (1.0 + std::log(diff_cumsum / delta_t));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_1D::dualEval(double point, double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  double a = (cumsum[t] - cumsum[s]) / (t - s); // Sbar_st
  double b = (cumsum[s] - cumsum[r]) / (s - r); // Sbar_rs

  ///
  point = point * muMax(a,b); /// Point rescaling
  ///

  return (costRecord[s] - minCost_t) / (t - s) + point * (costRecord[s] - costRecord[r]) / (s - r)
  + (1 - point) * 0.5 * (std::log((a - point * b) / (1 - point)) + 1);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_1D::dualMax(double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Variance_1D::muMax(double a, double b) const
{
  if (b != 0) return std::min(1., a/b);
  return 1.;
}

double Variance_1D::xMax(double a, double b) const
{
  if (a < b) return -a/(a-b);
  return std::numeric_limits<double>::infinity();
}

bool Variance_1D::isLeftBoundary(double a) const {return a == 0;}
double Variance_1D::Dstar_leftboundary() const {return std::numeric_limits<double>::infinity();}
double Variance_1D::Dstar_superLinearLimit() const {return 0;}

double Variance_1D::Dstar(double x) const
{
  return -0.5 * (std::log(x) + 1.0);
}


double Variance_1D::DstarPrime(double x) const
{
  return -0.5 /x;
}

double Variance_1D::DstarPrimeInv(double x) const
{
  return -0.5 /x;
}

double Variance_1D::DstarSecond(double x) const
{
  return 0.5/std::pow(x,2);
}

std::string Variance_1D::get_model() const { return "variance"; }
