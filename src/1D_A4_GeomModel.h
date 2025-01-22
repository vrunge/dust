#ifndef Geom_1D_H
#define Geom_1D_H

#include <Rcpp.h>

#include "1D_A_DUST.h"

using namespace Rcpp;

class Geom_1D : public DUST_1D {
public:
  Geom_1D(int dual_max, bool random_constraint, Nullable<int> nbLoops = Nullable<int>());
protected:
  double Cost(unsigned int t, unsigned int s) const override;
  double statistic(double& data) const override;

  double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const override;
  double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const override;

  double muMax(double a, double b) const override;

  double Dstar(double x) const override;
  double DstarPrime(double x) const override;
  double DstarSecond(double x) const override;

  std::string get_model() const override;
};

#endif
