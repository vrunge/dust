#include <cmath>

#include "MD_A3_ExpModel.h"

using namespace Rcpp;

Exp_MD::Exp_MD(int dual_max, bool random_constraint, Nullable<unsigned> nbLoops)
  : DUST_MD(dual_max, random_constraint, nbLoops) {}

double Exp_MD::Cost(const unsigned int& t, const unsigned int& s) const
{
  double diff;
  double res = 0;
  double inv_delta = pow(t - s, -1);
  for (unsigned int row = 0; row < d; row++)
  {
    diff = cumsum(row, t) - cumsum(row, s);
    if(diff <= 0){diff = 1e-100;} /// choice  1e-100 to avoid -Inf
    res += (1.0 + std::log(diff * inv_delta));
  }
  return res * (t - s);
}

double Exp_MD::statistic(const double& data) const
{return(data);}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Exp_MD::muMax(const double& a, const double& b) const
{
  double res = 1;
  if(b != 0){res = std::min(1.0, a/b);}
  return res;
}

std::array<double, 2> Exp_MD::muInterval(const arma::colvec& a, const arma::colvec& b, double& c, double& D) const
{
  std::array<double, 2> interval = {0, std::numeric_limits<double>::infinity() };

  for (unsigned int i = 0; i < a.n_elem; ++i)
  {
    if (a[i] > 0 && b[i] < 0){interval[1] = std::min(interval[1], -a[i]/b[i]);}
    else if (a[i] < 0 && b[i] > 0) {interval[0] = std::max(interval[0], -a[i]/b[i]);}
  }

  if (c > 0 && D < 0)
  {interval[1] = std::min(interval[1], -c / D);}
  else if (c < 0 && d > 0){interval[0] = std::max(interval[0], -c / D);}

  return(interval);
}

void Exp_MD::clipStepSizeModel(const double& m_elem, const arma::rowvec& constraint_means, const double& mu_sum, const arma::rowvec& direction, const double& direction_sum, double& max_stepsize) const
{
  // DIRECTION_SUM SHOULD BE POSITIVE HERE
  double dot_product = arma::dot(constraint_means, direction);
  if (dot_product > 0)
    return;

  double new_stepsize = -(1 + mu_sum) * m_elem / dot_product;
  if (new_stepsize < max_stepsize)
    max_stepsize = new_stepsize;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double Exp_MD::dual1D_Eval(double& point, const arma::colvec& a, const arma::colvec& b, double& c, double& d, double& e, double& f) const
{
  return(-std::numeric_limits<double>::infinity());
}

////////////////////////////////////////////////////////////////////////////////

double Exp_MD::dual1D_Max(double& argmax, arma::colvec& a, arma::colvec& b, double& c, double& D, double& e, double& f) const
{
  double Max = -std::numeric_limits<double>::infinity();

  return(Max);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double Exp_MD::Dstar(const double& x) const
{
  return (-std::log(x) - 1.0);
}


double Exp_MD::DstarPrime(const double& x) const
{
  return -1.0/x;
}

double Exp_MD::DstarSecond(const double& x) const
{
  return 1.0/std::pow(x,2);
}

std::string Exp_MD::get_model() const { return "exp"; }


