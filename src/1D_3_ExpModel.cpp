#include <Rcpp.h>
#include <cmath>

#include "1D_3_ExpModel.h"

using namespace Rcpp;

Exp_1D::Exp_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops)
  : DUST_1D(dualmax_algo, constr_index, nbLoops) {}

double Exp_1D::Cost(unsigned int t, unsigned int s) const
{
  double delta_t = t - s;
  double diff_cumsum = cumsum[t] - cumsum[s];
  if(diff_cumsum <= 0){diff_cumsum = 1e-100;} /// choice  1e-100 to avoid -Inf
  return delta_t * (1.0 + std::log(diff_cumsum / delta_t));
}

double Exp_1D::statistic(double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Exp_1D::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  double a = (cumsum[t] - cumsum[s]) / (t - s); // Sbar_st
  double b = (cumsum[s] - cumsum[r]) / (s - r); // Sbar_rs

  ///
  /// point in the right interval:
  if(b != 0){point = point * std::min(1.0, a/b);}
  ///
  ///

  return (costRecord[s] - minCost) / (t - s)
  + point * (costRecord[s] - costRecord[r]) / (s - r)
  + (1 - point) * (std::log((a - point * b) / (1 - point)) + 1);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Exp_1D::dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const
{
  return (-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Exp_1D::muMax(double a, double b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}

bool Exp_1D::isBoundary(double a) const
{
  return false;
}

double Exp_1D::Dstar(double x) const
{
  return (-std::log(x) - 1.0);
}


double Exp_1D::DstarPrime(double x) const
{
  return -1.0/x;
}

double Exp_1D::DstarPrimeInv(double x) const
{
  return -1.0/x;
}

double Exp_1D::DstarSecond(double x) const
{
  return 1.0/std::pow(x,2);
}

std::string Exp_1D::get_model() const
{
  return "exp";
}
