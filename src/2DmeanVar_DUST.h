#ifndef DUST_meanVar_H
#define DUST_meanVar_H

#include <Rcpp.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "1D_Indices.h"

using namespace Rcpp;

class DUST_meanVar
{
  public:
    DUST_meanVar(std::string dualmax_algo, std::string constr_index);

    virtual ~DUST_meanVar();

    void append_data(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());
    void update_partition();
    List get_partition();
    List get_info();
    List dust(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());

    ////////////////////////////////
    ////////////////////////////////
    ////////////////////////////////

  private:

    std::vector<double> cumsum;
    std::vector<double> cumsum2;
    std::vector<double> costRecord;

    double costEval(double point1, double point2, unsigned int t, unsigned int s);
    double costMin(unsigned int t, unsigned int s);

    double decisionEval1(double point, double minCost_t,
                         unsigned int t, unsigned int s, unsigned int r);
    double decisionEval2(double point1, double point2, double minCost_t,
                         unsigned int t, unsigned int s, unsigned int r);

    double dualEval(double point, double minCost,
                    unsigned int t, unsigned int s, unsigned int r);

    double dualMax(double minCost_t,
                   unsigned int t, unsigned int s, unsigned int r);

    double muMax(double a, double b, double a2, double b2);
    double xMax(double a, double b);

    bool isLeftBoundary(double a);
    double Dstar_leftboundary();
    double Dstar_superLinearLimit();

    double Dstar(double x);
    double DstarPrime(double x);
    double DstarPrimeInv(double x);
    double DstarSecond(double x);

    //////////// RANDOM NUMBER GENERATOR ////////////

    std::minstd_rand0 engine;  // Random number engine
    std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

    // --- // Test and Indices init // --- //
    void pruning_method();

    bool isOnePointOrLinear(double a, double b);
    bool specialCasePruning(double a,
                            double b,
                            double c,
                            double d,
                            double mu_max);

    // DUAL = -(c - mu * d)  - (1- mu) * Dstar((a - mu * b) / (1 - mu));

    // a = (cumsum[t] - cumsum[s]) / (t - s)
    // b = (cumsum[s] - cumsum[r]) / (s - r)
    // c = (Qt - Qs) / (t - s) =  (minCost_t - costRecord[s]) / (t - s);
    // d = (Qs - Qr) / (s - r)  = (costRecord[s] - costRecord[r]) / (s - r);

    ////////// MAX DUAL METHODS //////////
    ////////// MAX DUAL METHODS //////////
    ////////// MAX DUAL METHODS //////////
    // 0: DUSTr. random eval
    // 1: DUST. exact eval
    bool dualMaxAlgoRand(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxAlgoExact(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualOP(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxPelt(double minCost_t, unsigned int t, unsigned int s, unsigned int r);

    bool (DUST_meanVar::*current_test)(double minCost_t, unsigned int t,
                                                  unsigned int s,
                                                  unsigned int r);

    double penalty;

    std::string dualmax_algo;
    std::string constr_index;
    Indices_1D* indices;
    std::vector<int> nb_indices;
    unsigned int n; // number of observations

    ////////// Result processing //////////
    std::forward_list<unsigned int> backtrack_changepoints();
    std::vector<int> chptRecord;
};

#endif
