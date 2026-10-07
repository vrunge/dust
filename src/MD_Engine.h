/// MD DUST: dynamic programming (pruned and OP loops)

#ifndef MD_ENGINE_H
#define MD_ENGINE_H

#include <string>

#include "MD_Decision.h"
#include "MD_CoordinateDescent.h"
#include "MD_QuasiNewton.h"
#include "MD_RandomEval.h"
#include "MD_GaussExact.h"
#include "MD_Exact.h"

namespace dust_md {

class MDBase
{
public:
  MDBase(std::string model, std::string method, std::string backend,
         int constraints, int nb_iterations, double epsilon)
    : model_(std::move(model)), method_(std::move(method)),
      constraints_requested_(constraints), nb_iterations_(nb_iterations),
      epsilon_(epsilon)
  {
#ifdef HAVE_HIGHWAY
    highway_ = backend == "highway";
#else
    (void)backend;
#endif
  }
  virtual ~MDBase() = default;
  virtual void append_data(const NumericMatrix& data, Nullable<double> penalty) = 0;
  virtual void update_partition() = 0;
  virtual List get_partition() = 0;
  virtual List get_info() = 0;
  List dust(const NumericMatrix& data, Nullable<double> penalty)
  {
    append_data(data, penalty);
    update_partition();
    return get_partition();
  }

protected:
  std::string model_;
  std::string method_;
  int constraints_requested_;
  int nb_iterations_;
  double epsilon_;
  bool highway_ = false;
};

template <class Model>
class MDCommon : public MDBase
{
public:
  using MDBase::MDBase;

  void append_data(const NumericMatrix& data, Nullable<double> supplied) override
  {
    const size_t rows = data.nrow();
    const size_t cols = data.ncol();
    if (rows == 0) stop("data must have at least one row");
    if (dimension_ != 0 && rows != dimension_)
      stop("the number of rows cannot change after the first append");
    if (constraints_requested_ > static_cast<int>(rows))
      stop("constraints must be between 1 and the number of rows");
    if (cols == 0) return;
    for (double value : data) Math<Model>::validate(value);

    const double new_penalty = supplied.isNull() ?
      (dimension_ == 0 ? 2.0 * rows * std::log(static_cast<double>(cols)) : penalty_) :
      as<double>(supplied);
    if (!std::isfinite(new_penalty) || new_penalty < 0.0)
      stop("penalty must be a finite nonnegative number");
    if (dimension_ != 0 && new_penalty != penalty_)
      stop("penalty cannot change after the first nonempty append");

    if (dimension_ == 0)
    {
      dimension_ = rows;
      constraints_ = constraints_requested_ == 0 ? rows : constraints_requested_;
      penalty_ = new_penalty;
      sums_.assign(rows, 0.0);
      if constexpr (kVariance) sums_lo_.assign(rows, 0.0);
      costs_.push_back(-penalty_);
      previous_.push_back(0);
      active_.push_back(0);
#ifdef HAVE_HIGHWAY
      if (highway_)
      {
        highway_state_.model = dustib::model_id<Model>;
        highway_state_.sums.resize(rows);
        if constexpr (kVariance)
        {
          highway_state_.sums_lo.resize(rows);
          highway_state_.push_variance(0, costs_[0], sums_.data(), sums_lo_.data());
        }
        else highway_state_.push(0, costs_[0], sums_.data());
      }
#endif
    }
    sums_.reserve((length_ + cols + 1) * dimension_);
    if constexpr (kVariance)
    {
      // sums of y^2 + rounding errors (TwoSum)
      sums_lo_.reserve((length_ + cols + 1) * dimension_);
      for (size_t col = 0; col < cols; ++col)
        for (size_t row = 0; row < dimension_; ++row)
        {
          const size_t at = (length_ + col) * dimension_ + row;
          const double x = Model::statistic(data(row, col));
          const double a = sums_[at], sum = a + x, b = sum - a;
          sums_.push_back(sum);
          sums_lo_.push_back(sums_lo_[at] + ((a - (sum - b)) + (x - b)));
        }
    }
    else
      for (size_t col = 0; col < cols; ++col)
        for (size_t row = 0; row < dimension_; ++row)
          sums_.push_back(sums_[(length_ + col) * dimension_ + row] +
                          Model::statistic(data(row, col)));
    length_ += cols;
  }

  List get_partition() override
  {
    if (length_ == 0) stop("append data before requesting a partition");
    if (done_ != length_) stop("update_partition before requesting a partition");
    std::vector<unsigned int> changepoints{static_cast<unsigned int>(length_)};
    for (unsigned int at = previous_[length_]; at != 0; at = previous_[at])
      changepoints.push_back(at);
    std::reverse(changepoints.begin(), changepoints.end());
    std::vector<unsigned int> remaining(active_.rbegin(), active_.rend());
    return List::create(
      _["changepoints"] = changepoints,
      _["lastIndexSet"] = remaining,
      _["backend"] = highway_ ? "highway" : "scalar",
      _["nb"] = counts_,
      _["costQ"] = NumericVector(costs_.begin() + 1, costs_.end())
    );
  }

  List get_info() override
  {
    RObject epsilon_info = R_NilValue;
    if (epsilon_ >= 0.0) epsilon_info = wrap(epsilon_);
    return List::create(
      _["backend"] = highway_ ? "highway" : "scalar",
      _["data_length"] = length_,
      _["dimension"] = dimension_,
      _["current_penalty"] = penalty_,
      _["model"] = model_,
      _["pruning_algo"] = method_,
      _["constraints"] = constraints_,
      _["nbIterations"] = nb_iterations_,
      _["epsilon"] = epsilon_info
    );
  }

protected:
  std::pair<double, unsigned int> scan(size_t t,
                                       std::vector<double>* candidates) const
  {
#ifdef HAVE_HIGHWAY
    if (highway_)
    {
      if constexpr (kVariance)
      {
        const auto result = HWY_DYNAMIC_DISPATCH(md_hw::ScanVariance)(
          highway_state_, sums_.data() + t * dimension_,
          sums_lo_.data() + t * dimension_, static_cast<double>(t), candidates);
        return {result.first, active_[result.second]};
      }
      const auto result = HWY_DYNAMIC_DISPATCH(md_hw::Scan)(
        highway_state_, sums_.data() + t * dimension_,
        static_cast<double>(t), candidates);
      return {result.first, active_[result.second]};
    }
#endif
    if (candidates) candidates->resize(active_.size());
    double best = std::numeric_limits<double>::infinity();
    unsigned int argmin = 0;
    for (size_t i = 0; i < active_.size(); ++i)
    {
      const unsigned int s = active_[i];
      const double candidate = costs_[s] + cost(t, s);
      if (candidates) (*candidates)[i] = candidate;
      if (candidate < best) { best = candidate; argmin = s; }
    }
    return {best, argmin};
  }

  void retain_candidates(const std::vector<unsigned int>& survivors)
  {
#ifdef HAVE_HIGHWAY
    if (highway_)
    {
      if constexpr (kVariance) highway_state_.compact_variance(survivors);
      else highway_state_.compact(survivors);
    }
#else
    (void)survivors;
#endif
  }

  void push_candidate(size_t t)
  {
#ifdef HAVE_HIGHWAY
    if (highway_)
    {
      if constexpr (kVariance)
        highway_state_.push_variance(t, costs_[t], sums_.data() + t * dimension_,
                                     sums_lo_.data() + t * dimension_);
      else highway_state_.push(t, costs_[t], sums_.data() + t * dimension_);
    }
#else
    (void)t;
#endif
  }

  /// VARIANCE model: sums with rounding errors in sums_lo_
  static constexpr bool kVariance = dustib::model_id<Model> == 7;

  // sum of the statistic over (s, t] for dimension `row`
  double segment_sum(size_t t, size_t s, size_t row) const
  {
    const double plain = sums_[t * dimension_ + row] - sums_[s * dimension_ + row];
    if constexpr (kVariance)
      return plain + (sums_lo_[t * dimension_ + row] - sums_lo_[s * dimension_ + row]);
    else
      return plain;
  }

  double cost(size_t t, size_t s) const
  {
    const double span = static_cast<double>(t - s);
    double total = 0.0;
    for (size_t row = 0; row < dimension_; ++row)
    {
      double mean = segment_sum(t, s, row) / span;
      if constexpr (dustib::model_id<Model> == 4 || dustib::model_id<Model> == 5)
        mean = std::clamp(mean, 0.0, 1.0);
      total -= span * Math<Model>::conjugate(mean);
    }
    return total;
  }

  // decision function for s at time t, with `count` constraints
  const Decision<Model>& decision(size_t t, size_t s,
                                  const unsigned int* last, size_t count)
  {
    const unsigned int* earlier = last - count;
    Decision<Model>& test = test_;
    test.dimension = dimension_;
    test.constraints = count;
    test.highway = highway_;
    test.a.resize(dimension_);
    test.scratch.resize(dimension_);
    test.matrix.resize(dimension_ * count);
    test.u.resize(count);
    const double span = static_cast<double>(t - s);
    test.c = (costs_[t] - costs_[s]) / span;
    for (size_t row = 0; row < dimension_; ++row)
      test.a[row] = segment_sum(t, s, row) / span;
    for (size_t j = 0; j < count; ++j)
    {
      const size_t r = earlier[j];
      const double earlier_span = static_cast<double>(s - r);
      const double prior_cost = (costs_[s] - costs_[r]) / earlier_span;
      test.u[j] = test.c - prior_cost;
      for (size_t row = 0; row < dimension_; ++row)
      {
        const double prior_mean = segment_sum(s, r, row) / earlier_span;
        test.matrix[j * dimension_ + row] = test.a[row] - prior_mean;
      }
    }
    return test;
  }

  size_t dimension_ = 0;
  size_t constraints_ = 0;
  size_t length_ = 0;
  size_t done_ = 0;
  double penalty_ = 0.0;
  std::vector<double> sums_;
  std::vector<double> sums_lo_;   // variance model only
  std::vector<double> costs_;
  std::vector<unsigned int> previous_;
  std::vector<unsigned int> active_; // increasing order
  std::vector<int> counts_;
  Decision<Model> test_;                 // workspace
  std::vector<double> candidates_;
  std::vector<unsigned int> survivors_;
#ifdef HAVE_HIGHWAY
  md_hw::State highway_state_;
#endif
};

template <class Model>
class MDPruned : public MDCommon<Model>
{
public:
  using MDCommon<Model>::MDCommon;

  void update_partition() override
  {
    for (size_t t = this->done_ + 1; t <= this->length_; ++t)
    {
      std::vector<double>& candidates = this->candidates_;
      const auto best = this->scan(t, &candidates);
      this->costs_.push_back(best.first + this->penalty_);
      this->previous_.push_back(best.second);

      std::vector<unsigned int>& survivors = this->survivors_;
      survivors.clear();
      survivors.reserve(this->active_.size() + 1);
      for (size_t position = 0; position < this->active_.size(); ++position)
      {
        const unsigned int s = this->active_[position];
        const double candidate = candidates[position];
        const double scale = std::abs(this->costs_[s]) +
          std::abs(candidate - this->costs_[s]) + std::abs(this->costs_[t]);
        bool prune = dustib::positive(candidate - this->costs_[t], scale);
        // constraints = largest active indices smaller than s
        if (!prune && this->method_ != "PELT" && !survivors.empty())
        {
          const bool exact = this->method_ == "exact";
          size_t count = std::min(survivors.size(), this->constraints_);
          // exact, not gauss: 2 constraints at most
          if (exact && dustib::model_id<Model> != 0) count = std::min<size_t>(count, 2);
          const auto& test = this->decision(
            t, s, survivors.data() + survivors.size(), count);
          // no pruning certificate at x = 0 (not for gauss exact, k <= 2)
          const bool closed_form = exact && dustib::model_id<Model> == 0 && count <= 2;
          if (!closed_form && no_pruning_certificate(test, test.a.data()))
            prune = false;
          else if (exact)
            prune = exact_search(test);
          else if (this->method_ == "coordinateDescent")
            prune = coordinate_search(test, this->nb_iterations_,
                                      this->epsilon_);
          else if (this->method_ == "QN")
            prune = qn_search(test, this->nb_iterations_, this->epsilon_);
          else
            prune = random_search(test, this->nb_iterations_);
        }
        if (!prune)
        {
          survivors.push_back(s);
        }
      }
      this->counts_.push_back(static_cast<int>(survivors.size()));
      this->retain_candidates(survivors);
      survivors.push_back(static_cast<unsigned int>(t));
      this->active_.swap(survivors);
      this->push_candidate(t);
      this->done_ = t;
    }
  }
};

////////////////////////////////////////////////////////////////////////////////
/// OP: no pruning
template <class Model>
class MDOP : public MDCommon<Model>
{
public:
  using MDCommon<Model>::MDCommon;

  void update_partition() override
  {
    for (size_t t = this->done_ + 1; t <= this->length_; ++t)
    {
      const auto best = this->scan(t, nullptr);
      this->costs_.push_back(best.first + this->penalty_);
      this->previous_.push_back(best.second);
      this->counts_.push_back(static_cast<int>(this->active_.size()));
      this->active_.push_back(static_cast<unsigned int>(t));
      this->push_candidate(t);
      this->done_ = t;
    }
  }
};

} // namespace dust_md

#endif
