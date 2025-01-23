#include <cmath>

#include "MD_A5_BernModel.h"

using namespace Rcpp;

Bern_MD::Bern_MD(int dual_max, bool random_constraint, Nullable<int> nbLoops)
  : DUST_MD(dual_max, random_constraint, nbLoops) {}

double Bern_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double res = 0;
  double ratio;
  for (unsigned int row = 0; row < d; row++)
  {
    ratio = (cumsum(row, t) - cumsum(row, s)) / (t - s);
    if (ratio <= 0 | ratio >= 1)
      continue;
    res += - double(t - s) * (ratio * std::log(ratio) + (1 - ratio) * std::log(1 - ratio));
  }
  return res;
}



double Bern_MD::statistic(const double& data) const
{return(data);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Bern_MD::muMax(const double& a, const double& b) const
{
  double res = 1;
  if(b != 0 && b != 1){res = std::min(a/b, (1 - a)/(1 - b));}
  else{
    if(b == 0){res = 1 - a;}else{res = a;}
  }
  return res;
}

void Bern_MD::clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const
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

double Bern_MD::Dstar(const double& x) const
{
  return x*std::log(x) + (1.0 - x)*std::log(1.0 - x);
}


double Bern_MD::DstarPrime(const double& x) const
{
  return std::log(x) - std::log(1.0 - x);
}


double Bern_MD::DstarSecond(const double& x) const
{
  return 1.0/x + 1.0/(1.0 - x);
}

std::string Bern_MD::get_model() const { return "bern"; }







