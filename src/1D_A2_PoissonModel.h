#ifndef Poisson_1D_H
#define Poisson_1D_H

#include <Rcpp.h>

#include "1D_A_DUST.h"

using namespace Rcpp;

class Poisson_1D : public DUST_1D {
public:
  Poisson_1D(bool use_dual_max, bool random_constraint, Nullable<double> alpha = Nullable<double>());
protected:
  double Cost(unsigned int t, unsigned int s) const override;
  double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r) const override;
  double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r) const override;

};

#endif