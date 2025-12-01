#ifndef Gauss_1D_H
#define Gauss_1D_H

#include <Rcpp.h>

#include "1D_DUST.h"

using namespace Rcpp;

class Gauss_1D : public DUST_1D {
public:
  Gauss_1D(std::string dualmax_algo, std::string constr_index, Nullable<int> nbLoops = Nullable<int>());
protected:

  double statistic(double& data) const override;

  double costEval(double point, unsigned int t, unsigned int s) const override;
  double costMin(unsigned int t, unsigned int s) const override;

  double dualEval(double point, double minCost_t, unsigned int t, unsigned int s, unsigned int r) const override;
  double dualMax(double minCost_t, unsigned int t, unsigned int s, unsigned int r) const override;

  double muMax(double a, double b) const override;
  double xMax(double a, double b) const override;

  bool isLeftBoundary(double a) const override;
  double Dstar_leftboundary() const override;
  double Dstar_superLinearLimit() const override;


  double Dstar(double x) const override;
  double DstarPrime(double x) const override;
  double DstarPrimeInv(double x) const override;
  double DstarSecond(double x) const override;

  std::string get_model() const override;
};

#endif
