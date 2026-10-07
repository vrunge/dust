#ifndef DUST_1D_VARIANCE_H
#define DUST_1D_VARIANCE_H

/// VARIANCE model (statistic y^2): own scalar engines
///
/// With plain cumulative sums, S_t - S_s loses the small values y^2 when S_t
/// is large (e.g. y ~ 1e-7 after 1e6 data points) and we get false single
/// point segments. Here we also keep the rounding errors of the sums.
/// Otherwise same code as in 1D_DUST_Impl.h and 1D_OP_Impl.h.

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <type_traits>

#include "1D_DUST.h"
#include "1D_Indices.h"
#include "1D_DualMaxPolicies.h"
#include "1D_Validation.h"
#include "1D_8_VarianceModel.h"

using namespace Rcpp;

namespace dust_variance {

using Model = VariancePolicy;
constexpr int K = dustib::model_id<Model>;

/// cost of a segment (as VariancePolicy::costMin)
inline double segment_cost(double statistic, double length)
{
  if (statistic <= 0) statistic = 1e-100;
  return 0.5 * length * (1.0 + std::log(statistic / length));
}

/// pruning test with a, b, c, d (see 1D_DualMaxPolicies.h)
template <class Policy>
bool prunes(double a, double b, double c, double d)
{
  if constexpr (std::is_same<Policy, DualMax_DUST<Model>>::value)
    return DualMax_DUST<Model>::test_abcd(a, b, c, d);
  else if constexpr (std::is_same<Policy, DualMax_PELT<Model>>::value)
    return -Model::Dstar(a) - c > 0;
  else
    return dustib::test<K, 2>(a, b, c, d);   // DUSTib
}

////////////////////////////////////////////////////////////////////////////////
/// common part: sums in `cumsum`, rounding errors in `cumsum_lo_`
class VarianceBase : public DUST_1D
{
  public:
    explicit VarianceBase(std::string name) : DUST_1D(std::move(name)) {}

  protected:
    std::vector<double> cumsum_lo_;

    // S_t - S_s
    double statistic(unsigned int t, unsigned int s) const
    {
      return (cumsum[t] - cumsum[s]) + (cumsum_lo_[t] - cumsum_lo_[s]);
    }

    double cost(unsigned int t, unsigned int s) const
    {
      return segment_cost(statistic(t, s), static_cast<double>(t - s));
    }

    // returns true at the first append
    bool append_sums(const Rcpp::NumericVector& inData, Nullable<double> inPenalty)
    {
      const std::size_t m = static_cast<std::size_t>(inData.size());
      for (double y : inData) dustib::Math<K>::validate(y);
      const bool first_execution = (n == 0);
      const double checked_penalty = dust_validation::penalty(
        inPenalty, penalty, m, !first_execution);
      n += m;

      cumsum.reserve(n + 1);
      cumsum_lo_.reserve(n + 1);
      chptRecord.reserve(n + 1);
      costRecord.reserve(n + 1);

      if (first_execution)
      {
        penalty = checked_penalty;
        cumsum.push_back(0);
        cumsum_lo_.push_back(0);
        costRecord.push_back(-penalty);
        chptRecord.push_back(0);
      }
      for (std::size_t i = 0; i < m; ++i)
      {
        // TwoSum: s + error = a + x
        const double x = Model::statistic(inData[i]);
        const double a = cumsum.back();
        const double s = a + x;
        const double b = s - a;
        cumsum.push_back(s);
        cumsum_lo_.push_back(cumsum_lo_.back() + ((a - (s - b)) + (x - b)));
      }
      return first_execution;
    }

    std::string get_model() const override { return Model::get_model(); }
};

////////////////////////////////////////////////////////////////////////////////
/// DUST, DUSTib, PELT (as DUST_1D_T)
template <class Policy>
class Engine : public VarianceBase
{
  public:
    Engine() : VarianceBase(Policy::name()) {}

    void append_data(const Rcpp::NumericVector& inData, Nullable<double> inPenalty) override
    {
      if (inData.size() == 0) return;
      nb_indices.reserve(n + inData.size() + 1);
      if (append_sums(inData, inPenalty))
      {
        nb_indices.push_back(1);
        index_.add_first(0);
      }
    }

    void update_partition() override
    {
      if (n == 0) return;
      double lastCost;
      unsigned int nbt = nb_indices.back();
      for (unsigned t = index_.get_first() + 1; t <= n; t++)
      {
        // OP step
        index_.reset();
        double minCost_t = std::numeric_limits<double>::infinity();
        unsigned int argMin = 0;
        do
        {
          unsigned int s = index_.get_current();
          lastCost = costRecord[s] + cost(t, s); // without the penalty beta
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

        // pruning
        index_.reset_pruning();
        while (index_.is_not_the_last_pruning())
        {
          const unsigned int s = index_.get_current();
          const unsigned int r = index_.get_constraint();
          const double a = statistic(t, s) / (t - s);
          const double b = statistic(s, r) / (s - r);
          const double c = (minCost_t - costRecord[s]) / (t - s);
          const double d = (costRecord[s] - costRecord[r]) / (s - r);
          if (prunes<Policy>(a, b, c, d))
          {
            index_.prune_current();
            nbt--;
          }
          else
          {
            index_.next_pruning();
          }
        }
        // last index: PELT test
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

  private:
    Indices_1D index_;
};

////////////////////////////////////////////////////////////////////////////////
/// OP (as DUST_1D_OP_T)
class EngineOP : public VarianceBase
{
  public:
    EngineOP() : VarianceBase("OP") {}

    void append_data(const Rcpp::NumericVector& inData, Nullable<double> inPenalty) override
    {
      if (inData.size() == 0) return;
      append_sums(inData, inPenalty);
    }

    void update_partition() override
    {
      for (unsigned t = last_done + 1; t <= n; t++)
      {
        double minCost_t = std::numeric_limits<double>::infinity();
        unsigned int argMin = 0;
        for (unsigned int s = 0; s < t; s++)
        {
          double c = costRecord[s] + cost(t, s);
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
      std::forward_list<unsigned int> lastIndexSet;
      for (unsigned int i = 0; i <= n; i++) lastIndexSet.push_front(i);
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

  private:
    unsigned int last_done = 0;
};

////////////////////////////////////////////////////////////////////////////////
inline DUST_1D* make_engine(const std::string& method)
{
  if (method == "DUST")   return new Engine<DualMax_DUST<Model>>();
  if (method == "DUSTib") return new Engine<DualMax_DUSTib<Model>>();
  if (method == "PELT")   return new Engine<DualMax_PELT<Model>>();
  if (method == "OP")     return new EngineOP();
  stop("Unrecognized method \"" + method + "\". Valid values: \"DUST\", \"DUSTib\", \"PELT\", \"OP\".");
}

} // namespace dust_variance

#endif
