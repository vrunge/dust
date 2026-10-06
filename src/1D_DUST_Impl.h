#ifndef DUST_1D_IMPL_H
#define DUST_1D_IMPL_H

#include <Rcpp.h>
#include <cmath>

#include "1D_DUST.h"
#include "1D_Indices.h"
#include "1D_DualMaxPolicies.h"
#include "1D_Validation.h"

using namespace Rcpp;

// Scalar segmentation engine for a given model and pruning rule.

template <class Model, class DualMaxPolicy>
class DUST_1D_T : public DUST_1D
{
  public:
    DUST_1D_T() : DUST_1D(DualMaxPolicy::name()) {}

    // Add observations and update the cumulative sufficient statistics.

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
      nb_indices.reserve(n + 1);

      if (first_execution)
      {
        penalty = checked_penalty;

        cumsum.push_back(0);
        costRecord.push_back(-penalty);
        chptRecord.push_back(0);
        nb_indices.push_back(1);

        index_.add_first(0);
      }

      for (std::size_t i = 0; i < m; ++i)
        cumsum.push_back(cumsum.back() + Model::statistic(inData[i]));
    }

    void update_partition() override
    {
      if (n == 0) return;
      double lastCost;
      unsigned int nbt = nb_indices.back();
      for (unsigned t = index_.get_first() + 1; t <= n; t++)
      {
        // Optimal partitioning step.
        index_.reset();
        double minCost_t = std::numeric_limits<double>::infinity();
        unsigned int argMin = 0;
        do
        {
          unsigned int s = index_.get_current();
          lastCost = costRecord[s] + Model::costMin(cumsum, t, s); // without the penalty beta
          if (lastCost < minCost_t)
          {
            minCost_t = lastCost;
            argMin = s;
          }
          index_.next();
        }
        while(index_.is_not_the_last());
        minCost_t += penalty;
        costRecord.push_back(minCost_t);
        chptRecord.push_back(argMin);

        // Apply the pruning rule.
        index_.reset_pruning();

        while (index_.is_not_the_last_pruning())
        {
          if (DualMaxPolicy::test(cumsum, costRecord, minCost_t, t, index_.get_current(), index_.get_constraint()))
          {
            index_.prune_current();
            nbt--;
          }
          else
          {
            index_.next_pruning();
          }
        }
        // The last index gives the PELT test (zero multiplier).
        if (lastCost > minCost_t)
        {
          index_.prune_last();
          nbt--;
        }

        index_.add_first(t);
        nb_indices.push_back(nbt);
        nbt++;
      }
    }

    List get_partition() override
    {
      if (n == 0) stop("append data before requesting a partition");
      if (chptRecord.size() != static_cast<std::size_t>(n) + 1)
        stop("update_partition before requesting a partition");
      std::forward_list<unsigned int> chpts = backtrack_changepoints();

      return List::create(
        _["changepoints"] = chpts,
        _["lastIndexSet"] = index_.get_list(),
        _["backend"] = "scalar",
        _["nb"] = Rcpp::IntegerVector(nb_indices.begin() + 1, nb_indices.end()),
        _["costQ"] = Rcpp::NumericVector(costRecord.begin() + 1, costRecord.end())
      );
    }

  protected:
    std::string get_model() const override { return Model::get_model(); }

  private:
    Indices_1D index_;
};

#endif
