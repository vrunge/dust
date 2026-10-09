#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include <hwy/contrib/thread_pool/thread_pool.h>

#include "coordinate_search.h"
#include "decision_1d.h"
#include "detector.h"
#include "exact_search.h"
#include "qn_search.h"

namespace dust {

Detector::Detector(Method method, Cost cost, double penalty, Options options)
  : method_(method), cost_(cost), penalty_(penalty), options_(std::move(options)),
    scan_(scan_function(cost, method == Method::OP || method == Method::FrontPELT))
{
  const Options& o = options_;
  if (!std::isfinite(penalty) || penalty < 0.0)
    throw std::invalid_argument("penalty must be a finite nonnegative number");
  if (o.constraints < 0) throw std::invalid_argument("constraints must be positive or NULL");
  if (o.iterations < 1) throw std::invalid_argument("nbIterations must be a positive integer");
  if (o.epsilon != -1.0 && (!std::isfinite(o.epsilon) || o.epsilon < 0.0))
    throw std::invalid_argument("epsilon must be finite, nonnegative, or NULL");
  if (o.epsilon >= 0.0 && o.solver == Solver::RandomEval)
    throw std::invalid_argument("epsilon is not available for randomEval");
  if (o.num_threads < 1) throw std::invalid_argument("the number of threads must be a positive integer");
  if (cost == Cost::GaussianMeanVariance && method == Method::DUST &&
      !(o.solver == Solver::Exact && (o.constraints == 1 || o.constraints == 2)))
    throw std::invalid_argument("GaussianMeanVariance: DUST with the Exact solver and 1 or 2 constraints");
  if (!options_.uniform)
    options_.uniform = [generator = std::mt19937_64()]() mutable {
      return std::uniform_real_distribution<double>()(generator);
    };
}

////////////////////////////////////////////////////////////////////////////////
/// one time step of the dynamic programming for the cost K
template <Cost K>
class Detector::Step
{
public:
  Step(Detector& s, hwy::ThreadPool* pool) : s_(s), pool_(pool) {}

  void operator()(size_t t)
  {
    const size_t arg = scan(t);
    const double qt = val_[arg] + s_.penalty_;
    s_.costs_.push_back(qt);
    s_.previous_.push_back(static_cast<int>(s_.candidates_.position(arg)));
    prune(t);
    s_.nb_.push_back(static_cast<int>(s_.candidates_.size()));
    s_.candidates_.push(t, qt, s_.sums_);
  }

private:
  /// val_ = Q_s + C(s, t) for the candidates, returns the argmin. A run of
  /// the pool costs about 1 microsecond per thread, the threads share the
  /// candidates from kPerThread candidates per thread (4 tasks per thread for
  /// the work stealing between cores of different speeds)
  size_t scan(size_t t)
  {
    constexpr size_t kPerThread = 4096;
    const Candidates& c = s_.candidates_;
    const size_t k = c.size(), threads = pool_ ? pool_->NumWorkers() : 1;
    val_.resize(k);
    if (threads == 1 || k < kPerThread * threads) return s_.scan_(c, s_.sums_, t, 0, k, val_.data());
    const size_t tasks = 4 * threads;
    argmins_.resize(tasks);
    pool_->Run(0, tasks, [&](uint64_t task, size_t) {
      argmins_[task] = s_.scan_(c, s_.sums_, t, task * k / tasks, (task + 1) * k / tasks, val_.data());
    });
    size_t best = argmins_[0];
    for (size_t arg : argmins_)
      if (val_[arg] < val_[best]) best = arg;
    return best;
  }

  /// dual test after the PELT test
  enum class Dual { None, ClosedForm, ClosedFormIB, MeanVariance, Search };

  void prune(size_t t)
  {
    switch (s_.method_)
    {
      case Method::OP: return;
      case Method::FrontPELT: return prune_front(t);
      case Method::PELT: return prune<Dual::None>(t);
      case Method::DUST: break;
    }
    if constexpr (K == Cost::GaussianMeanVariance) prune<Dual::MeanVariance>(t);
    else
    {
      const Solver solver = s_.options_.solver;
      if (s_.n_features_ == 1 && solver == Solver::Exact) prune<Dual::ClosedForm>(t);
      else if (solver == Solver::ExactIB) prune<Dual::ClosedFormIB>(t);
      else prune<Dual::Search>(t);
    }
  }

  /// PELT test on the smallest candidates only
  void prune_front(size_t t)
  {
    size_t i = 0;
    while (i < s_.candidates_.size() && pelt_prunes(t, i)) ++i;
    s_.candidates_.drop_front(i);
  }

  /// kept candidates moved to [0, kept): the constraints of the dual tests
  template <Dual M>
  void prune(size_t t)
  {
    Candidates& c = s_.candidates_;
    const size_t k = c.size();
    // one feature: sums of the statistics (0, and 1 for GaussianMeanVariance)
    const double *pos = c.positions(), *cost = c.costs(), *hi = c.sums_hi(0), *lo = c.sums_lo(0);
    const double *hi_t = s_.sums_.hi_at(t), *lo_t = s_.sums_.lo_at(t);
    const double qt = s_.costs_.data()[t], tt = static_cast<double>(t);
    size_t kept = 0;
    for (size_t i = 0; i < k; ++i)
    {
      if (pelt_prunes(t, i)) continue;
      if constexpr (M == Dual::ClosedForm || M == Dual::ClosedFormIB)
        if (kept > 0)
        {
          const size_t r = kept - 1;
          const double inv_span = 1.0 / (tt - pos[i]), inv_before = 1.0 / (pos[i] - pos[r]);
          const double a = ((hi_t[0] - hi[i]) + (lo_t[0] - lo[i])) * inv_span;
          const double b = ((hi[i] - hi[r]) + (lo[i] - lo[r])) * inv_before;
          const double cc = (qt - cost[i]) * inv_span, d = (cost[i] - cost[r]) * inv_before;
          if (M == Dual::ClosedForm ? closed_form_prunes<K>(a, b, cc, d) : closed_form_ib_prunes<K>(a, b, cc, d))
            continue;
        }
      // a segment of one point has an infinite cost
      if constexpr (M == Dual::MeanVariance)
        if (kept > 0 && tt - pos[i] >= 2)
        {
          const double *hi2 = c.sums_hi(1), *lo2 = c.sums_lo(1);
          const double inv = 1.0 / (tt - pos[i]);
          const Moments a{((hi_t[0] - hi[i]) + (lo_t[0] - lo[i])) * inv,
                          ((hi_t[1] - hi2[i]) + (lo_t[1] - lo2[i])) * inv, (qt - cost[i]) * inv};
          const auto before = [&](size_t r) {
            const double inv = 1.0 / (pos[i] - pos[r]);
            return Moments{((hi[i] - hi[r]) + (lo[i] - lo[r])) * inv,
                           ((hi2[i] - hi2[r]) + (lo2[i] - lo2[r])) * inv, (cost[i] - cost[r]) * inv};
          };
          const Moments b = before(kept - 1);
          if (s_.options_.constraints == 1 ? mean_variance_prunes(a, b)
                                           : mean_variance_prunes(a, b, kept >= 2 ? before(kept - 2) : b))
            continue;
        }
      if constexpr (M == Dual::Search)
        if (kept > 0 && search_prunes(t, i, kept)) continue;
      if (kept != i) c.move(i, kept);
      ++kept;
    }
    c.resize(kept);
  }

  /// Q_s + C(s, t) > Q_t (beyond rounding errors)
  bool pelt_prunes(size_t t, size_t i) const
  {
    const double value = val_.data()[i], qt = s_.costs_.data()[t];
    if (!(value > qt)) return false;
    const double qs = s_.candidates_.cost_at(i);
    return positive(value - qt, std::abs(qs) + std::abs(value - qs) + std::abs(qt));
  }

  /// maximum of the decision function, the last kept candidates as constraints
  bool search_prunes(size_t t, size_t i, size_t kept)
  {
    const Detector& s = s_;
    const Candidates& cand = s.candidates_;
    const Options& o = s.options_;
    const size_t at = cand.position(i), dim = s.n_features_;
    const double qs = cand.cost_at(i), span = static_cast<double>(t - at);
    // exact, not gauss: 2 constraints at most
    const bool exact = o.solver == Solver::Exact;
    size_t count = std::min(kept, o.constraints == 0 ? dim : static_cast<size_t>(o.constraints));
    if (exact && K != Cost::GaussianMean) count = std::min<size_t>(count, 2);
    Decision<K>& test = decision_;
    test.resize(dim, count);
    // prefix sums of the features at one time are contiguous
    const PrefixSums& sums = s.sums_;
    const double *hi_s = sums.hi_at(at), *lo_s = sums.lo_at(at), *hi_t = sums.hi_at(t), *lo_t = sums.lo_at(t);
    test.c = (s.costs_[t] - qs) / span;
    for (size_t row = 0; row < dim; ++row)
      test.a[row] = ((hi_t[row] - hi_s[row]) + (lo_t[row] - lo_s[row])) / span;
    for (size_t j = 0; j < count; ++j)
    {
      const size_t r = cand.position(kept - count + j);
      const double *hi_r = sums.hi_at(r), *lo_r = sums.lo_at(r), before = static_cast<double>(at - r);
      test.u[j] = test.c - (qs - cand.cost_at(kept - count + j)) / before;
      for (size_t row = 0; row < dim; ++row)
        test.matrix[j * dim + row] = test.a[row] - ((hi_s[row] - hi_r[row]) + (lo_s[row] - lo_r[row])) / before;
    }

    // no pruning certificate at x = 0 (not for gauss exact, k <= 2: closed form)
    if (!(exact && K == Cost::GaussianMean && count <= 2) && no_pruning_certificate(test, test.a.data()))
      return false;
    switch (o.solver)
    {
      case Solver::CoordinateDescent: return coordinate_search(test, o.iterations, o.epsilon);
      case Solver::QN: return qn_search(test, o.iterations, o.epsilon);
      case Solver::RandomEval: return random_search(test, o.iterations, o.uniform);
      default: return exact_search(test);
    }
  }

  Detector& s_;
  hwy::ThreadPool* pool_;
  std::vector<double> val_;      // Q_s + C(s, t) for the candidates
  std::vector<size_t> argmins_;  // of the scan tasks
  Decision<K> decision_;
};

Detector& Detector::fit(const double* X, size_t n_samples, size_t n_features)
{
  n_samples_ = n_features_ = 0;
  costs_.clear();
  previous_.clear();
  nb_.clear();
  chgpts_.clear();
  return partial_fit(X, n_samples, n_features);
}

Detector& Detector::partial_fit(const double* X, size_t n_samples, size_t n_features)
{
  if (n_features == 0) throw std::invalid_argument("data must have at least one feature");
  if (n_features_ != 0 && n_features != n_features_)
    throw std::invalid_argument("the number of features cannot change after the first fit");
  if (cost_ != Cost::GaussianMeanVariance && options_.constraints > static_cast<int>(n_features))
    throw std::invalid_argument("constraints must be between 1 and the number of features");
  if (method_ == Method::DUST && n_features != 1 &&
      (cost_ == Cost::GaussianMeanVariance || options_.solver == Solver::ExactIB))
    throw std::invalid_argument("this solver and cost are for one feature");
  if (n_samples == 0) return *this;

  with_cost(cost_, [&](auto k) {
    constexpr Cost K = decltype(k)::value;
    using F = Family<K>;
    for (size_t i = 0; i < n_samples * n_features; ++i) F::validate(X[i]);
    const size_t series = n_features * F::statistics;
    if (n_samples_ == 0)
    {
      sums_.reset(series);
      costs_.push_back(-penalty_);
      previous_.push_back(0);
      candidates_.reset(series);
      candidates_.push(0, -penalty_, sums_);
    }
    std::vector<double> statistics(series);
    for (size_t t = 0; t < n_samples; ++t)
    {
      for (size_t j = 0; j < n_features; ++j) F::statistic(X[t * n_features + j], &statistics[j * F::statistics]);
      sums_.add(statistics.data());
    }
    n_samples_ += n_samples;
    n_features_ = n_features;

    // threads for the scans during this call (the caller is one of them)
    hwy::AlignedUniquePtr<hwy::ThreadPool> pool;
    if (options_.num_threads > 1) pool = hwy::MakeUniqueAligned<hwy::ThreadPool>(options_.num_threads - 1);
    Step<K> step(*this, pool.get());
    for (size_t t = costs_.size(); t <= n_samples_; ++t) step(t);
  });
  if (!std::isfinite(costs_.back())) throw std::domain_error("no finite segmentation is available for these data");

  chgpts_.clear();
  for (int t = previous_[n_samples_]; t != 0; t = previous_[t]) chgpts_.push_back(t);
  std::reverse(chgpts_.begin(), chgpts_.end());
  return *this;
}

std::vector<int> Detector::candidates() const
{
  std::vector<int> left;
  for (size_t i = fitted(candidates_).size(); i-- > 0;) left.push_back(static_cast<int>(candidates_.position(i)));
  return left;
}

} // namespace dust
