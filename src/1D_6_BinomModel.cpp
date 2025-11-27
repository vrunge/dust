#include <Rcpp.h>
#include <cmath>

#include "1D_6_BinomModel.h"

using namespace Rcpp;

Binom_1D::Binom_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops)
  : DUST_1D(dualmax_algo, constr_index, nbLoops) {}

double Binom_1D::Cost(unsigned int t, unsigned int s) const
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


double Binom_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double a = (cumsum[t] - cumsum[s]) / (t - s); // Sbar_st
  double b = (cumsum[s] - cumsum[r]) / (s - r); // Sbar_rs

  ///
  /// point in the right interval:
  if(b != 0 && b != 1){point = point * std::min(a/b, (1 - a)/(1 - b));}
  else{
    if(b == 0){point = point * (1 - a);}else{point = point * a;}
  }
  ///
  ///
  double R = (a - point * b) / (1 - point);

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  - (1 - point) * (R * std::log(R) + (1 - R) * std::log(1 - R));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Binom_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
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

bool Binom_1D::isBoundary(double a) const
{
  if(a == 0 || a == 1){return true;}
  return false;
}


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
