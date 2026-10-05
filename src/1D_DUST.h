#ifndef DUST_1D_H
#define DUST_1D_H

#include <Rcpp.h>
#include <forward_list>

using namespace Rcpp;

// --------------------------------------------------------------------- //
// DUST_1D: R-facing abstract base class.
//
// Holds everything that does not depend on the cost model or the pruning
// algorithm (bookkeeping, the methods exposed by the Rcpp module), and
// declares as pure virtual the few entry points that do: append_data,
// update_partition, get_partition, get_model. Those are implemented
// once, in the class template DUST_1D_T<Model, DualMaxPolicy> (see
// 1D_DUST_Impl.h): each cost model is a header-only policy struct of
// static functions, each pruning algorithm likewise (1D_DualMaxPolicies.h),
// both held or called by concrete type -- no vtable indirection, and (as
// of DualMaxPolicy becoming a template parameter rather than a runtime
// member-function-pointer choice) no indirect-call boundary either,
// anywhere in the OP-step scan or the pruning loop, which is where
// virtually all the running time goes. The active-index set (Indices_1D,
// see 1D_Indices.h) is likewise a plain, non-virtual, concretely-typed
// member of DUST_1D_T: there used to be a choice of index-selection
// strategy templated alongside Model, but the alternative (random
// constraint choice) measured strictly worse -- slower per pruning test
// and no stronger a pruning rule -- and was removed, so there is nothing
// left to select between.
//
// append_data/dust take a NumericVector by const reference rather than a
// std::vector<double>: Rcpp binds a NumericVector parameter directly to
// the R-owned SEXP, with no copy, whereas binding a std::vector<double>
// forces Rcpp to allocate and copy the whole input first.
// --------------------------------------------------------------------- //

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

    ////////////////////////////////
    ////////////////////////////////
    ////////////////////////////////

  protected:
    std::vector<double> cumsum;
    std::vector<double> costRecord;

    // Called once per get_info() call (negligible frequency): no need to
    // template this one.
    virtual std::string get_model() const = 0;

    std::string dualmax_algo;
    std::vector<int> nb_indices;
    unsigned int n; // number of observations
    double penalty;

    ////////////////////////////////
    ////////////////////////////////
    ////////////////////////////////

    ////////// Result processing //////////
    std::forward_list<unsigned int> backtrack_changepoints();
    std::vector<int> chptRecord;
};

#endif
