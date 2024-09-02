#ifndef DUST_meanVar_H
#define DUST_meanVar_H

#include <Rcpp.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "1D_B_Indices.h"

using namespace Rcpp;


class DUST_meanVar
{

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

public:
  DUST_meanVar(bool use_dual_max,
               bool random_constraint,
               Nullable<double> alpha = Nullable<double>(),
               Nullable<int> nbLoops = Nullable<int>());

  ~DUST_meanVar();

  // --- // Setup // --- //
  // fit is accessible by user
  void init(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());

  // --- // Main computation // --- //
  void compute(std::vector<double>& inData);

  // --- // Result retrieval // --- //
  // get_partition is accessible by user
  List get_partition();

  // --- // Wrapper method for quick use of the class // --- //
  // quick is accessible by user
  List quick(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

  private:
  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio

  std::vector<double> cumsum;
  std::vector<double> cumsum2;
  std::vector<double> costRecord;
  int nb_Loops; // number of loops in optimization step (For dual max)

  double Cost(unsigned int t, unsigned int s);
  double dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r);
  double dualMax(double minCost, unsigned int t, unsigned int s, unsigned int r);

  double Dstar(double x) const;
  double DstarPrime(double x) const;
  double DstarSecond(double x) const;

  //////////// RANDOM NUMBER GENERATOR ////////////

  std::minstd_rand0 engine;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

  // --- // Test and Indices init // --- //
  void init_method();

  // --- // Test handling // --- //
  double exact_test(double minCost, unsigned int t, unsigned int s, unsigned int r);
  double random_test(double minCost, unsigned int t, unsigned int s, unsigned int r);

  double (DUST_meanVar::*current_test)(double minCost, unsigned int t, unsigned int s, unsigned int r);

  // --- // Result processing // --- //
  std::forward_list<unsigned int> backtrack_changepoints();

  // --- // Private fields // --- //
  bool use_dual_max;
  bool random_constraint;
  double alpha;

  Indices* indices;
  std::forward_list<unsigned int> nb_indices;

  unsigned int n; // number of observations
  double penalty;

  std::vector<int> changepointRecord;

};

#endif
