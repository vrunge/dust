#include <Rcpp.h>
#include <algorithm> // for std::min
#include <cmath>

#include <random>

#include "1D_2_PoissonModel.h"

using namespace Rcpp;

Poisson_1D::Poisson_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops)
  : DUST_1D(dualmax_algo, constr_index, nbLoops) {}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_1D::costEval(double point, unsigned int t, unsigned int s) const
{
  return (std::exp(point) - point*(cumsum[t] - cumsum[s])/(1.0*(t-s)));
}

double Poisson_1D::costMin(unsigned int t, unsigned int s) const
{
  double res = 0.0;
  double diff = cumsum[t] - cumsum[s];
  if(diff > 0.0) /// should be != 0, we use > for robustness
  {
    res = diff * (1.0 - std::log(diff / (t - s)));
  }
  return res;
}

double Poisson_1D::statistic(double& data) const
{return(data);}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_1D::dualEval(double point, double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  double a = (cumsum[t] - cumsum[s]) / (t - s); // Sbar_st
  double b = (cumsum[s] - cumsum[r]) / (s - r); // Sbar_rs

  ///
  point = point * muMax(a,b); /// Point rescaling
  ///

  return (costRecord[s] - minCost_t) / (t - s) + point * (costRecord[s] - costRecord[r]) / (s - r)
    - (a - point * b) * (std::log((a - point * b) / (1 - point)) - 1);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_1D::dualMax(double minCost_t, unsigned int t, unsigned int s, unsigned int r) const
{
  ///
  /// -Dstar(a - (b-a)*x) - (c - (d-c)*x)
  ///
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost_t - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double mu_max = muMax(a, b);

  double ratio = -(c-d)/(a - b);
  double val = DstarPrimeInv(ratio);

  double x = 1.0/(a - b) * (val - a);
  double x_max = mu_max/(1 - mu_max);
  if(x < 0){x = 0;}
  //Rcout << " xxxxx: " <<  x << std::endl;
  if(x > x_max)
  {
    //Rcout << " MUmax: " << mu_max << " xmax: " <<  x_max << " x: " <<  x << std::endl;
    x = 0;
  }

  //Rcout << -Dstar(a + (a-b)*x) + d + (c + d)*x;
  //Rcout <<  " /a-b/ " << a - b;
  //Rcout << " /EVAL/ " << a + (a-b)*x;
  //Rcout << " /firstX/ " << 1.0/(a - b) * (val - a);
  //Rcout << " /VAL/ " << val;
  //Rcout << " /RATIO/ " << ratio;
  //Rcout << " //a: " << a << " b: " << b << " x: " << x << " mu_max: "<< mu_max << " mu_max: "<< mu_max << std::endl;
  //x = 0;
  return -Dstar(a + (a-b)*x) - (c + (c-d)*x);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Poisson_1D::muMax(double a, double b) const
{
  if (b != 0) return std::min(1., a/b);
  return 1.;
}

double Poisson_1D::xMax(double a, double b) const
{
  if (a < b) return -a/(a-b);
  else  return std::numeric_limits<double>::infinity();
}


bool Poisson_1D::isLeftBoundary(double a) const {return a == 0;}
double Poisson_1D::Dstar_leftboundary() const {return 0;}
double Poisson_1D::Dstar_superLinearLimit() const {return std::numeric_limits<double>::infinity();};

double Poisson_1D::Dstar(double x) const
{
  return (x * (std::log(x) - 1.0));
}


double Poisson_1D::DstarPrime(double x) const
{
  return std::log(x);
}

double Poisson_1D::DstarPrimeInv(double x) const
{
  return std::exp(x);
}

double Poisson_1D::DstarSecond(double x) const
{
  return pow(x, -1);
}

std::string Poisson_1D::get_model() const { return "poisson"; }
