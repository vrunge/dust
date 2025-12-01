#include <Rcpp.h>
#include <cmath>

#include "1D_4_GeomModel.h"

using namespace Rcpp;

Geom_1D::Geom_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops)
  : DUST_1D(dualmax_algo, constr_index, nbLoops) {}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Geom_1D::costEval(double point, unsigned int t, unsigned int s) const
{
  return -(t-s)*std::log(std::exp(-point) - 1) - point*(cumsum[t] - cumsum[s]);
}

double Geom_1D::costMin(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s])/(t - s);
  if(m > 1)
  {res = (t - s) * std::log(m - 1) - (cumsum[t] - cumsum[s]) * std::log((m - 1) / m);}
  return res;
}


double Geom_1D::statistic(double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Geom_1D::dualEval(double point, double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  double a = (cumsum[t] - cumsum[s]) / (t - s); // Sbar_st
  double b = (cumsum[s] - cumsum[r]) / (s - r); // Sbar_rs

  ///
  point = point * muMax(a,b); /// Point rescaling
  ///

  double R = (a - point * b) / (1 - point);

  return (costRecord[s] - minCost_t) / (t - s)  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + (1 - point) * (R * std::log(R) - (R - 1) * std::log(R - 1));
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Geom_1D::dualMax(double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Geom_1D::muMax(double a, double b) const
{
  if (b != 1) return std::min(1.0, (a-1)/(b-1));
  return 1.;
}


double Geom_1D::xMax(double a, double b) const
{
  if (a < b) return -(a-1)/(a-b);
  else  return std::numeric_limits<double>::infinity();
}


bool Geom_1D::isLeftBoundary(double a) const {return a == 1;}
double Geom_1D::Dstar_leftboundary() const {return 0;}
double Geom_1D::Dstar_superLinearLimit() const {return 0;}

double Geom_1D::Dstar(double x) const
{
  return (x - 1.0)*std::log(x - 1.0) - x*std::log(x);
}


double Geom_1D::DstarPrime(double x) const
{
  return std::log(x - 1.0) - std::log(x);
}

double Geom_1D::DstarPrimeInv(double x) const
{
  return 1/(1 - std::exp(x));
}

double Geom_1D::DstarSecond(double x) const
{
  return 1.0/(x - 1.0) - 1.0/x;
}

std::string Geom_1D::get_model() const { return "geom"; }
