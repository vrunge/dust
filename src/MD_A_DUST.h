#ifndef DUST_MD_H
#define DUST_MD_H

#include <RcppArmadillo.h>
#include <dqrng.h>
#include <xoshiro.h>

#include <fstream>

#include "MD_B_Indices.h"

#include "utils.h"

using namespace Rcpp;


class DUST_MD
{

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

public:
  DUST_MD(int dual_max,
          bool random_constraint,
          Nullable<double> alpha = Nullable<double>(),
          Nullable<int> nbLoops = Nullable<int>());

  virtual ~DUST_MD();

  // --- // Setup // --- //
  // fit is accessible by user
  void init(const arma::dmat& inData, Nullable<double> inPenalty = Nullable<double>(), Nullable<unsigned int> inNbMax = Nullable<unsigned int>());

  // --- // Main computation // --- //
  void compute(const arma::dmat& inData);

  // --- // Result retrieval // --- //
  // get_partition is accessible by user
  List get_partition();

  // --- // Wrapper method for quick use of the class // --- //
  // quick is accessible by user
  List quick(const arma::dmat& inData, Nullable<double> inPenalty = Nullable<double>(), Nullable<unsigned int> inNbMax = Nullable<unsigned int>());

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

protected:
  unsigned int n; // number of observations
  unsigned int d;
  unsigned int nb_max;

  const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
  const double m1 = 0.01;  // Armijo
  arma::dmat cumsum;
  std::vector<double> costRecord;
  int nb_Loops; // number of loops in optimization step (For dual max)

  virtual double Cost(const unsigned int& t, const unsigned int& s) const = 0;
  virtual double statistic(const double& value) const = 0;

  virtual double muMax(const double& a, const double& b) const = 0;

  virtual double dualMax(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r) = 0;
  virtual double dualEval(std::vector<unsigned int> point, const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r) = 0;

  virtual double Dstar(const double& x) const = 0;
  virtual double DstarPrime(const double& x) const = 0;
  virtual double DstarSecond(const double& x) const = 0;

  //////////// RANDOM NUMBER GENERATOR ////////////

  dqrng::xoshiro256plus engine;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////

private:
  // --- // Test and Indices init // --- //
  void init_method();

  // --- // MAX DUAL METHODS // --- //
  // --- //   // --- //   // --- //   // --- //
  bool dualMaxAlgo0(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r);
  bool dualMaxAlgo1(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r);
  // bool dualMaxAlgo3(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r);
  bool dualMaxAlgo4(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r);
  bool dualMaxAlgo5(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r);

  bool (DUST_MD::*current_test)(const double& minCost, const unsigned int& t, const unsigned int& s, std::vector<unsigned int> r);

  // --- // Result processing // --- //
  std::forward_list<unsigned int> backtrack_changepoints();

  // --- // Private fields // --- //
  int dual_max;
  bool random_constraint;
  double alpha;

  Indices_MD* indices;
  std::vector<int> nb_indices;
  double penalty;

  std::vector<int> changepointRecord;

  arma::rowvec mu;
  arma::rowvec mu_max;
  arma::rowvec inv_max;
  arma::rowvec grad;

  arma::rowvec linearTerm;

  arma::colvec m_value; // m(mu)
  arma::colvec objectiveMean;
  arma::dmat constraintMean;
  arma::colvec nonLinearGrad;

  arma::dmat Identity;
  arma::dmat inverseHessian; // for quasi newton optimizing

};

#endif
