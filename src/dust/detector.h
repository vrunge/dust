/// Optimal partitioning with pruning (Runge, Truong and Querne 2025, arXiv:2507.02467)
/// X: n_samples x n_features, row major

#ifndef DUST_DETECTOR_H
#define DUST_DETECTOR_H

#include <stdexcept>
#include <vector>

#include "costs.h"
#include "random_search.h"
#include "scan.h"

namespace dust {

enum class Method { OP, PELT, FrontPELT, DUST };   // FrontPELT: PELT test on the smallest candidates only

/// maximization of the dual (DUST)
enum class Solver { Exact, ExactIB, CoordinateDescent, QN, RandomEval };

struct Options
{
  Solver solver = Solver::Exact;
  int constraints = 1;   // 0: one per feature
  int iterations = 1;
  double epsilon = -1;   // < 0: fixed number of iterations
  int num_threads = 1;
  Uniform uniform;       // default: std::mt19937_64
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

  const std::vector<int>& chgpts() const { return fitted(chgpts_); }        // 0 < tau < n_samples
  std::vector<double> min_costs() const;   // Q_t, t = 0..n_samples, costs: -2 log-likelihood
  const std::vector<int>& n_candidates() const { return fitted(nb_); }
  std::vector<int> candidates() const;
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
  PrefixSums sums_;
  Candidates candidates_;
  std::vector<double> costs_;
  std::vector<int> previous_;   // last change point before t
  std::vector<int> nb_;
  std::vector<int> chgpts_;
};

} // namespace dust

#endif
