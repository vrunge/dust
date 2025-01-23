#include <cmath>

#include "MD_A4_GeomModel.h"

using namespace Rcpp;

Geom_MD::Geom_MD(int dual_max, bool random_constraint, Nullable<int> nbLoops)
  : DUST_MD(dual_max, random_constraint, nbLoops) {}

double Geom_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double res = 0;
  double ratio;
  double diff;
  double delta = double(t - s);
  double inv_delta = pow(delta, -1);
  for (unsigned int row = 0; row < d; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    ratio = diff * inv_delta;
    if (ratio <= 1)
      continue;
    res += delta * std::log(ratio - 1) - diff * std::log((ratio - 1) / ratio);
  }
return res;
}


double Geom_MD::statistic(const double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Geom_MD::dualMax(arma::colvec& a, arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  return 0.;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Geom_MD::muMax(const double& a, const double& b) const
{
  double res = 1;
  if(b != 1){res = std::min(1.0, (a-1)/(b-1));}
  return res;
}

void Geom_MD::clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const
{
  double dot_product = arma::dot(constraint_means, direction);
  double sum_diff = dot_product - direction_sum;

  if (dot_product < 0)
  {
    double new_stepsize = -(1 + mu_sum) * m_elem / dot_product;
    if (new_stepsize < max_stepsize)
      max_stepsize = new_stepsize;
  }

  if (sum_diff > 0)
  {
    double new_stepsize = (1 + mu_sum) * (1 - m_elem) / sum_diff;
    if (new_stepsize < max_stepsize)
      max_stepsize = new_stepsize;
  }
}

double Geom_MD::Dstar(const double& x) const
{
  return (x - 1.0)*std::log(x - 1.0) - x*std::log(x);
}


double Geom_MD::DstarPrime(const double& x) const
{
  return std::log(x - 1.0) - std::log(x);
}

double Geom_MD::DstarSecond(const double& x) const
{
  return 1.0/(x - 1.0) - 1.0/x;
}

std::string Geom_MD::get_model() const { return "geom"; }
