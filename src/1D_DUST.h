#ifndef DUST_1D_H
#define DUST_1D_H

#include <Rcpp.h>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

#include "1D_Indices.h"

using namespace Rcpp;

class DUST_1D
{
  public:
    DUST_1D(std::string dualmax_algo,
            std::string constr_index,
            Nullable<int> nbLoops = Nullable<int>());

    virtual ~DUST_1D();

    void append_data(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());
    void update_partition();
    List get_partition();
    List get_info();
    List dust(std::vector<double>& inData, Nullable<double> inPenalty = Nullable<double>());

    ////////////////////////////////
    ////////////////////////////////
    ////////////////////////////////

  protected:
    const double phi = (1 + sqrt(5)) / 2;  // Golden ratio
    const double m1 = 0.01;  // Armijo
    std::vector<double> cumsum;
    std::vector<double> costRecord;
    int nb_Loops; // number of loops in optimization step (For dual max)

    virtual double statistic(double& data) const = 0;

    virtual double costEval(double point,
                            unsigned int t,
                            unsigned int s) const = 0;

    virtual double costMin(unsigned int t,
                           unsigned int s) const = 0;

    virtual double dualEval(double point,
                            double minCost_t,
                            unsigned int t,
                            unsigned int s,
                            unsigned int r) const = 0;
    virtual double dualMax(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r) const = 0;


    virtual double muMax(double a, double b) const = 0;
    virtual double xMax(double a, double b) const = 0;

    virtual bool isLeftBoundary(double a) const = 0;
    virtual double Dstar_leftboundary() const = 0;
    virtual double Dstar_superLinearLimit() const = 0;

    virtual double Dstar(double x) const = 0;
    virtual double DstarPrime(double x) const = 0;
    virtual double DstarPrimeInv(double x) const = 0;
    virtual double DstarSecond(double x) const = 0;

    virtual std::string get_model() const = 0;

    //////////// RANDOM NUMBER GENERATOR ////////////

    std::minstd_rand0 engine;  // Random number engine
    std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

    ////////////////////////////////
    ////////////////////////////////
    ////////////////////////////////

  private:
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
    // 2: DUSTgs. golden-section search
    // 3: DUSTbs. binary search. At each step, we evaluate the tangent line to the current point at its max to stop the search at early step (when possible)
    // 4: DUSTqn.  Quasi-Newton
    // 5: PELT
    // 6: OP
    // 7: DUST. exact eval
    bool dualMaxAlgo0(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxAlgo1(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxAlgo2(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxAlgo3(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxAlgo4(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxAlgo5(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxAlgo6(double minCost_t, unsigned int t, unsigned int s, unsigned int r);
    bool dualMaxAlgo7(double minCost_t, unsigned int t, unsigned int s, unsigned int r);

    bool (DUST_1D::*current_test)(double minCost_t, unsigned int t,
                                                  unsigned int s,
                                                  unsigned int r);

    double penalty;

    //int nb00 = 0;
    //int nb0 = 0;
    //int nb1 = 0;
    //int nb2 = 0;
    //int nb00T = 0;
    //int nb0T = 0;
    //int nb1T = 0;
    //int nb2T = 0;

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
