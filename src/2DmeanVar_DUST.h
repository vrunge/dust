#ifndef DUST_meanVar_H
#define DUST_meanVar_H

#include <Rcpp.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "2D_Indices.h"

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

    std::minstd_rand0 engine;  // Random number engine
    std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

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
                         unsigned int t, unsigned int s, unsigned int r1, unsigned int r2);


    double dualEval(double point, double minCost,
                    unsigned int t, unsigned int s, unsigned int r);

    double dualMax(double minCost_t,
                   unsigned int t, unsigned int s, unsigned int r);

    double muMax(double a, double b, double a2, double b2);
    double xMax(double a, double b);

    // --- // Test and Indices init // --- //
    void pruning_method();


    // DUAL = -(c - mu * d)  - (1- mu) * Dstar((a - mu * b) / (1 - mu));

    // a = (cumsum[t] - cumsum[s]) / (t - s)
    // b = (cumsum[s] - cumsum[r]) / (s - r)
    // c = (Qt - Qs) / (t - s) =  (minCost_t - costRecord[s]) / (t - s);
    // d = (Qs - Qr) / (s - r)  = (costRecord[s] - costRecord[r]) / (s - r);

    ////////// MAX DUAL METHODS //////////
    ////////// MAX DUAL METHODS //////////
    ////////// MAX DUAL METHODS //////////
    bool dualMaxAlgoRand(double minCost_t, unsigned int t, unsigned int s, unsigned int r1, unsigned int r2);
    bool dualOP(double minCost_t, unsigned int t, unsigned int s, unsigned int r1, unsigned int r2);
    bool dualMaxPelt(double minCost_t, unsigned int t, unsigned int s, unsigned int r1, unsigned int r2);
    bool decisionTest1(double minCost_t, unsigned int t, unsigned int s, unsigned int r1, unsigned int r2);
    bool decisionTest2(double minCost_t, unsigned int t, unsigned int s, unsigned int r1, unsigned int r2);

    bool (DUST_meanVar::*current_test)(double minCost_t, unsigned int t,
                                                  unsigned int s,
                                                  unsigned int r1,
                                                  unsigned int r2);

    double penalty;

    std::string dualmax_algo;
    std::string constr_index;
    Indices_2D* indices;
    std::vector<int> nb_indices;
    unsigned int n; // number of observations

    ////////// Result processing //////////
    std::forward_list<unsigned int> backtrack_changepoints();
    std::vector<int> chptRecord;
};

#endif
