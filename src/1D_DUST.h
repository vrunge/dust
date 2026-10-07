#ifndef DUST_1D_H
#define DUST_1D_H

#include <Rcpp.h>
#include <forward_list>

using namespace Rcpp;

/// scalar engines: common class (see DUST_1D_T)

class DUST_1D
{
  public:
    DUST_1D(std::string dualmax_algo);

    virtual ~DUST_1D();

    virtual void append_data(const Rcpp::NumericVector& inData, Nullable<double> inPenalty = Nullable<double>()) = 0;
    virtual void update_partition() = 0;
    virtual List get_partition() = 0;
    List get_info();
    List dust(const Rcpp::NumericVector& inData, Nullable<double> inPenalty = Nullable<double>());

  protected:
    std::vector<double> cumsum;
    std::vector<double> costRecord;

    virtual std::string get_model() const = 0;

    std::string dualmax_algo;
    std::vector<int> nb_indices;
    unsigned int n; // number of observations
    double penalty;

    std::forward_list<unsigned int> backtrack_changepoints();
    std::vector<int> chptRecord;
};

#endif
