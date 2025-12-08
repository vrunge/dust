#include <Rcpp.h>
#include <cmath>

#include "1D_6_BinomModel.h"

using namespace Rcpp;

Binom_1D::Binom_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops)
  : DUST_1D(dualmax_algo, constr_index, nbLoops) {}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Binom_1D::costEval(double point, unsigned int t, unsigned int s) const
{
  const double dt   = static_cast<double>(t) - static_cast<double>(s);
  const double mean = (cumsum[t] - cumsum[s]) / dt;
  double val = 0;
  if (point <= 0.0) { val = std::log1p(std::exp(point));}
  else { val =  point + std::log1p(std::exp(-point));}
  // Original: log(exp(point) + 1) - point * (cumsum[t] - cumsum[s]) / (1.0*(t-s))
  return val - point * mean;
}


double Binom_1D::costMin(unsigned int t, unsigned int s) const
{
  double res = 0;
  double m = (cumsum[t] - cumsum[s]) / (t - s);
  if(m > 0 && m < 1)
  {res = - double(t - s) * (m * std::log(m) + (1 - m) *  std::log(1 - m));}
  return res;
}

double Binom_1D::statistic(double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Binom_1D::dualEval(double point, double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
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

double Binom_1D::dualMax(double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Binom_1D::muMax(double a, double b) const
{
  double res = 1;
  if(b != 0 && b != 1){res = std::min(a/b, (1 - a)/(1 - b));}
  else{
    if(b == 0){res = 1 - a;}else{res = a;}
  }
  return res;
}

double Binom_1D::xMax(double a, double b) const
{
//  if(std::abs(a-b) < 1e-10){return ;}
  if (a < b) return -a/(a-b);
  else  return (1-a)/(a-b);
}


bool Binom_1D::isLeftBoundary(double a) const {return a == 0;}
double Binom_1D::Dstar_leftboundary() const {return 0;}
double Binom_1D::Dstar_superLinearLimit() const {return 0;}

double Binom_1D::Dstar(double x) const
{
  return x*std::log(x) + (1.0 - x)*std::log(1.0 - x);
}


double Binom_1D::DstarPrime(double x) const
{
  return std::log(x) - std::log(1.0 - x);
}

double Binom_1D::DstarPrimeInv(double x) const
{
  return std::exp(x) / (1 + std::exp(x));
}


double Binom_1D::DstarSecond(double x) const
{
  return 1.0/x + 1.0/(1.0 - x);
}

std::string Binom_1D::get_model() const { return "binom"; }
