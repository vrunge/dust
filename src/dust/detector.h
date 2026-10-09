/// Change point detection: optimal partitioning with pruning
/// (Runge, Truong and Querne 2025, arXiv:2507.02467), used as in the Python
/// package changepoints:
///
///   dust::Detector detector(dust::Method::DUST, dust::Cost::GaussianMean, penalty);
///   detector.fit(X, n_samples, n_features);
///   detector.chgpts();   // change points
///
/// X: n_samples x n_features, row major (X[t * n_features + j]); the features
/// are independent series with the same cost. partial_fit adds samples to the
/// fitted ones (online use).
///
/// Q_t = min over the candidates s of Q_s + C(s, t) + penalty, Q_0 = -penalty.
/// After each step, the method removes candidates that cannot be the last
/// change point at a later time:
///   OP          no pruning
///   PELT        Q_s + C(s, t) > Q_t
///   FrontPELT   PELT test on the smallest candidates only: the candidates stay
///               consecutive (faster scan, with threads), more of them are kept
///   DUST        PELT test, then the dual test: s is pruned if the dual is
///               positive somewhere, the constraints are the largest candidates
///               kept below s. Options::solver looks for that point.
/// Cost::GaussianMeanVariance: DUST with one feature and 1 or 2 constraints.
/// Options::num_threads: threads sharing the scans of the candidates (Highway
/// thread pool, for the duration of a fit).

#ifndef DUST_DETECTOR_H
#define DUST_DETECTOR_H

#include <stdexcept>
#include <vector>

#include "costs.h"
#include "random_search.h"
#include "scan.h"

namespace dust {

enum class Method { OP, PELT, FrontPELT, DUST };

/// maximization of the dual of DUST
enum class Solver
{
  Exact,               // exact maximum: closed form for one feature (decision_1d.h),
                       // exact_search.h otherwise (at most 2 constraints but for GaussianMean)
  ExactIB,             // one feature: closed form with explicit domain checks
  CoordinateDescent,   // coordinate ascent (coordinate_search.h)
  QN,                  // quasi-Newton (qn_search.h)
  RandomEval           // random evaluations (random_search.h)
};

struct Options
{
  Solver solver = Solver::Exact;   // DUST
  int constraints = 1;   // DUST: number of constraints (0: one per feature)
  int iterations = 1;    // CoordinateDescent, QN, RandomEval: number of iterations
  double epsilon = -1;   // CoordinateDescent, QN: precision instead (if >= 0)
  int num_threads = 1;   // scan of the candidates
  Uniform uniform;       // RandomEval: uniform random numbers (default: std::mt19937_64)
};

class Detector
{
public:
  Detector(Method method, Cost cost, double penalty, Options options = {});

  Detector& fit(const double* X, size_t n_samples, size_t n_features);
  Detector& partial_fit(const double* X, size_t n_samples, size_t n_features);
  const std::vector<int>& predict() const { return chgpts(); }
  std::vector<int> fit_predict(const double* X, size_t n_samples, size_t n_features)
  {
    return fit(X, n_samples, n_features).predict();
  }

  Method method() const { return method_; }
  Cost cost() const { return cost_; }
  double penalty() const { return penalty_; }
  const Options& options() const { return options_; }

  /// fitted attributes
  const std::vector<int>& chgpts() const { return fitted(chgpts_); }        // change points, 0 < tau < n_samples
  const std::vector<double>& min_costs() const { return fitted(costs_); }   // Q_t, t = 0..n_samples
  const std::vector<int>& n_candidates() const { return fitted(nb_); }      // candidates kept at t = 1..n_samples
  std::vector<int> candidates() const;                                      // candidates left (decreasing)
  size_t n_samples() const { return n_samples_; }
  size_t n_features() const { return n_features_; }

private:
  template <Cost K> class Step;
  template <class T>
  const T& fitted(const T& attribute) const
  {
    if (n_samples_ == 0) throw std::logic_error("fit the detector first");
    return attribute;
  }

  Method method_;
  Cost cost_;
  double penalty_;
  Options options_;
  ScanFunction scan_;
  size_t n_samples_ = 0, n_features_ = 0;
  PrefixSums sums_;              // statistics of the features
  Candidates candidates_;
  std::vector<double> costs_;    // Q_t
  std::vector<int> previous_;    // argmin: last change point before t
  std::vector<int> nb_;          // candidates kept at t
  std::vector<int> chgpts_;
};

} // namespace dust

#endif
