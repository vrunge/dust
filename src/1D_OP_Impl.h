#ifndef DUST_1D_OP_IMPL_H
#define DUST_1D_OP_IMPL_H

#include <Rcpp.h>
#include <cmath>
#include <limits>

#include "1D_DUST.h"
#include "1D_DUSTib.h"
#include "1D_Validation.h"

using namespace Rcpp;

// Optimal Partitioning scans every previous endpoint. It has its own
// recursion because the general DUST engine also checks its oldest active
// candidate for pruning.

template <class Model>
class DUST_1D_OP_T : public DUST_1D
{
  public:
    DUST_1D_OP_T() : DUST_1D("OP") {}

    void append_data(const Rcpp::NumericVector& inData, Nullable<double> inPenalty) override
    {
      const std::size_t m = static_cast<std::size_t>(inData.size());
      if (m == 0) return;
      for (double y : inData) dustib::Math<dustib::model_id<Model>>::validate(y);
      bool first_execution = (n == 0);
      const double checked_penalty = dust_validation::penalty(
        inPenalty, penalty, m, !first_execution);
      n += m;

      cumsum.reserve(n + 1);
      chptRecord.reserve(n + 1);
      costRecord.reserve(n + 1);

      if (first_execution)
      {
        penalty = checked_penalty;

        cumsum.push_back(0);
        costRecord.push_back(-penalty);
        chptRecord.push_back(0);
      }

      for (std::size_t i = 0; i < m; ++i)
        cumsum.push_back(cumsum.back() + Model::statistic(inData[i]));
    }

    void update_partition() override
    {
      for (unsigned t = last_done + 1; t <= n; t++)
      {
        double minCost_t = std::numeric_limits<double>::infinity();
        unsigned int argMin = 0;
        for (unsigned int s = 0; s < t; s++)
        {
          double c = costRecord[s] + Model::costMin(cumsum, t, s);
          if (c < minCost_t)
          {
            minCost_t = c;
            argMin = s;
          }
        }
        minCost_t += penalty;
        costRecord.push_back(minCost_t);
        chptRecord.push_back(argMin);
      }
      last_done = n;
    }

    List get_partition() override
    {
      if (n == 0) stop("append data before requesting a partition");
      if (chptRecord.size() != static_cast<std::size_t>(n) + 1)
        stop("update_partition before requesting a partition");
      std::forward_list<unsigned int> chpts = backtrack_changepoints();

      // Nothing is ever pruned: the full active set is {0, 1, ..., n},
      // in the same largest-to-smallest order Indices_1D::get_list() uses
      // (see 1D_Indices.h) -- built directly rather than via that class,
      // since there is no index_ member here at all.
      std::forward_list<unsigned int> lastIndexSet;
      for (unsigned int i = 0; i <= n; i++) lastIndexSet.push_front(i);

      // nb is trivially t for every t (no index 0 excluded, matching
      // DUST_1D_T's convention: nb_indices[0]=1 is never reported, and
      // nb_indices[t] for the t-th step counts {0,...,t-1}, i.e. t
      // entries, once the newly-added t itself is excluded).
      std::vector<int> nb(n);
      for (unsigned int t = 1; t <= n; t++) nb[t - 1] = static_cast<int>(t);

      return List::create(
        _["changepoints"] = chpts,
        _["lastIndexSet"] = lastIndexSet,
        _["backend"] = "scalar",
        _["nb"] = nb,
        _["costQ"] = Rcpp::NumericVector(costRecord.begin() + 1, costRecord.end())
      );
    }

  protected:
    std::string get_model() const override { return Model::get_model(); }

  private:
    unsigned int last_done = 0;
};

#endif
