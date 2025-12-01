#include <Rcpp.h>
#include <cmath>

#include "1D_5_BernModel.h"

using namespace Rcpp;

Bern_1D::Bern_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops)
  : DUST_1D(dualmax_algo, constr_index, nbLoops) {}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Bern_1D::costEval(double point, unsigned int t, unsigned int s) const
{
  return (t-s)*std::log(std::exp(point) + 1) - point*(cumsum[t] - cumsum[s]);
}

double Bern_1D::costMin(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s]) / (t - s);
  if(m > 0 && m < 1)
  {res = - double(t - s) * (m * std::log(m) + (1 - m) *  std::log(1 - m));}
  return res;
}

double Bern_1D::statistic(double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Bern_1D::dualEval(double point, double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  double a = (cumsum[t] - cumsum[s]) / (t - s); // Sbar_st
  double b = (cumsum[s] - cumsum[r]) / (s - r); // Sbar_rs

  ///
  point = point * muMax(a,b); /// Point rescaling
  ///

  double R = (a - point * b) / (1 - point);

  return (costRecord[s] - minCost_t) / (t - s) + point * (costRecord[s] - costRecord[r]) / (s - r)
  - (1 - point) * (R * std::log(R) + (1 - R) * std::log(1 - R));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Bern_1D::dualMax(double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Bern_1D::muMax(double a, double b) const
{
  double res = 1;
  if(b != 0 && b != 1){res = std::min(a/b, (1 - a)/(1 - b));}
  else{
    if(b == 0){res = 1 - a;}else{res = a;}
  }
  return res;
}

double Bern_1D::xMax(double a, double b) const
{
  if (a < b) return -a/(a-b);
  else  return (1-a)/(a-b);
}


bool Bern_1D::isLeftBoundary(double a) const {return a == 0;}
double Bern_1D::Dstar_leftboundary() const {return 0;}
double Bern_1D::Dstar_superLinearLimit() const {return 0;}

double Bern_1D::Dstar(double x) const
{
  return x*std::log(x) + (1.0 - x)*std::log(1.0 - x);
}

double Bern_1D::DstarPrime(double x) const
{
  return std::log(x) - std::log(1.0 - x);
}

double Bern_1D::DstarPrimeInv(double x) const
{
  return std::exp(x) / (1 + std::exp(x));
}

double Bern_1D::DstarSecond(double x) const
{
  return 1.0/x + 1.0/(1.0 - x);
}

std::string Bern_1D::get_model() const { return "bern"; }
