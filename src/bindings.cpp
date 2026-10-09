#include <Rcpp.h>
#include <optional>

#include "dust/detector.h"
#include "dust/preprocessing.h"

using namespace Rcpp;

namespace {

dust::Cost parse_model(const std::string& model)
{
  static const std::pair<const char*, dust::Cost> costs[] = {
    {"gauss", dust::Cost::GaussianMean}, {"poisson", dust::Cost::Poisson},
    {"exp", dust::Cost::Exponential}, {"geom", dust::Cost::Geometric},
    {"bern", dust::Cost::Bernoulli}, {"binom", dust::Cost::Binomial},
    {"negbin", dust::Cost::NegativeBinomial}, {"variance", dust::Cost::GaussianVariance},
    {"meanVar", dust::Cost::GaussianMeanVariance}};
  for (const auto& c : costs)
    if (model == c.first) return c.second;
  stop("Unrecognized model \"" + model + "\".");
}

/// meanVar "1D" and "2D": DUST with 1 or 2 constraints
std::pair<dust::Method, dust::Solver> parse_method(const std::string& method)
{
  using M = dust::Method;
  using S = dust::Solver;
  static const std::pair<const char*, std::pair<M, S>> methods[] = {
    {"OP", {M::OP, S::Exact}}, {"PELT", {M::PELT, S::Exact}}, {"PELTpar", {M::FrontPELT, S::Exact}},
    {"DUST", {M::DUST, S::Exact}}, {"DUSTib", {M::DUST, S::ExactIB}}, {"exact", {M::DUST, S::Exact}},
    {"coordinateDescent", {M::DUST, S::CoordinateDescent}}, {"QN", {M::DUST, S::QN}},
    {"randomEval", {M::DUST, S::RandomEval}}, {"1D", {M::DUST, S::Exact}}, {"2D", {M::DUST, S::Exact}}};
  for (const auto& m : methods)
    if (method == m.first) return m.second;
  stop("Unrecognized method \"" + method + "\".");
}

/// built at the first nonempty append: the default penalty depends on the first data
class RDetector
{
public:
  RDetector(std::string model, std::string method, int constraints, int iterations,
            double epsilon, int threads)
    : model_(std::move(model)), method_(std::move(method)), cost_(parse_model(model_)),
      parsed_(parse_method(method_))
  {
    options_ = {parsed_.second, method_ == "2D" ? 2 : constraints, iterations, epsilon, threads,
                [] { return unif_rand(); }};
    dust::Detector(parsed_.first, cost_, 0.0, options_);
  }

  void append_data(NumericVector data, Nullable<double> penalty)
  {
    size_t d = 1, n = data.size();
    if (data.hasAttribute("dim"))
    {
      const IntegerVector dim = data.attr("dim");
      d = dim[0];
      n = dim[1];
    }
    if (d == 0) stop("data must have at least one row");
    if (d_ != 0 && d != d_) stop("the number of rows cannot change after the first append");
    if (n == 0) return;
    if (cost_ == dust::Cost::GaussianMeanVariance && d != 1)
      stop("meanVar requires one feature");
    const double dof = cost_ == dust::Cost::GaussianMeanVariance ? 2.0 : static_cast<double>(d);
    if (!detector_)
    {
      const double value = penalty.isNull() ? 2.0 * dof * std::log(static_cast<double>(n)) : as<double>(penalty);
      detector_.emplace(parsed_.first, cost_, value, options_);
      d_ = d;
    }
    else if (!penalty.isNull() && as<double>(penalty) != detector_->penalty())
      stop("penalty cannot change after the first nonempty append");
    pending_.insert(pending_.end(), data.begin(), data.end());
    length_ += n;
  }

  void update_partition()
  {
    if (pending_.empty()) return;
    RNGScope scope;   // randomEval
    detector_->partial_fit(pending_.data(), pending_.size() / d_, d_);
    pending_.clear();
  }

  List get_partition() const
  {
    if (!detector_) stop("append data before requesting a partition");
    if (!pending_.empty()) stop("update_partition before requesting a partition");
    std::vector<int> changepoints = detector_->chgpts();
    changepoints.push_back(static_cast<int>(detector_->n_samples()));
    const std::vector<double> costs = detector_->min_costs();
    return List::create(
      _["changepoints"] = changepoints,
      _["lastIndexSet"] = detector_->candidates(),
      _["nb"] = detector_->n_candidates(),
      _["costQ"] = std::vector<double>(costs.begin() + 1, costs.end()));
  }

  List get_info() const
  {
    size_t constraints = options_.constraints == 0 ? d_ : static_cast<size_t>(options_.constraints);
    if (parsed_.first == dust::Method::DUST && parsed_.second == dust::Solver::Exact &&
        cost_ != dust::Cost::GaussianMean && cost_ != dust::Cost::GaussianMeanVariance)
      constraints = std::min<size_t>(constraints, 2);
    return List::create(
      _["data_length"] = length_,
      _["dimension"] = d_,
      _["current_penalty"] = detector_ ? detector_->penalty() : 0.0,
      _["model"] = model_,
      _["pruning_algo"] = method_,
      _["constraints"] = constraints,
      _["nbIterations"] = options_.iterations,
      _["epsilon"] = options_.epsilon >= 0.0 ? wrap(options_.epsilon) : R_NilValue,
      _["threads"] = options_.num_threads);
  }

  List dust(NumericVector data, Nullable<double> penalty)
  {
    append_data(data, penalty);
    update_partition();
    return get_partition();
  }

private:
  std::string model_, method_;
  dust::Cost cost_;
  std::pair<dust::Method, dust::Solver> parsed_;
  dust::Options options_;
  std::optional<dust::Detector> detector_;
  size_t d_ = 0, length_ = 0;
  std::vector<double> pending_;   // appended, not analysed yet
};

} // namespace

RCPP_MODULE(DUSTMODULE)
{
  class_<RDetector>("Detector")
    .constructor<std::string, std::string, int, int, double, int>()
    .method("append_data", &RDetector::append_data)
    .method("update_partition", &RDetector::update_partition)
    .method("get_partition", &RDetector::get_partition)
    .method("get_info", &RDetector::get_info)
    .method("dust", &RDetector::dust);
}

//' sdDiff
//'
//' @description Estimation of the noise standard deviation from the differences of the data (robust to the changes in mean)
//'
//' @param y a numeric vector
//' @param method \code{"HALL"} (default), \code{"MAD"} or \code{"SD"}
//' \itemize{
//'   \item \code{"HALL"}: Hall estimator (weighted differences)
//'   \item \code{"MAD"}: MAD of the differences
//'   \item \code{"SD"}: standard deviation of the differences
//' }
//' @return the estimated standard deviation
//'
//' @examples
//' y <- dataGenerator_1D(chpts = c(300, 600), parameters = c(0, 1), sdNoise = 2, type = "gauss")
//' sdDiff(y)
//' sdDiff(y, "MAD")
//' sdDiff(y, "SD")
//'
//' @export
// [[Rcpp::export]]
double sdDiff(std::vector<double> y, std::string method = "HALL")
{
  return dust::sd_diff(y, method);
}

//' data_normalization_1D
//'
//' @name data_normalization_1D
//'
//' @description Model-specific transformations of univariate data
//' \itemize{
//'   \item \code{"gauss"}: division by \code{sdDiff(y)}
//'   \item \code{"poisson"}, \code{"exp"}: division by the mean
//'   \item \code{"binom"}, \code{"negbin"}: division by \code{size}
//'   \item \code{"variance"}: the mean is removed
//'   \item \code{"geom"}, \code{"bern"}: no change
//' }
//'
//' @param y a numeric vector
//' @param type the model: \code{"gauss"} (default), \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}, \code{"variance"}
//' @param size number of trials (binom) or number of successes (negbin). Required for these two models.
//' @return the normalized data
//' @details For Poisson likelihood segmentation, use the original counts with the default penalty.
//' If the counts are divided by their mean, divide the penalty by the same mean to preserve the objective.
//' For binomial and negative binomial data, divide the penalty by \code{size} as well.
//' Exponential rescaling only adds a segmentation-independent constant to the cost.
//' @examples
//' y <- dataGenerator_1D(chpts = c(300, 600), parameters = c(0, 1), sdNoise = 2, type = "gauss")
//' sdDiff(data_normalization_1D(y))
//'
//' y <- dataGenerator_1D(chpts = c(300, 600), parameters = c(2, 8), type = "poisson")
//' mean(data_normalization_1D(y, type = "poisson"))
//'
//' y <- dataGenerator_1D(chpts = c(300, 600), parameters = c(0.4, 0.7), nbTrials = 5, type = "binom")
//' data_normalization_1D(y, type = "binom", size = 5)[1:10]
//'
//' @export
// [[Rcpp::export]]
std::vector<double> data_normalization_1D(std::vector<double> y, std::string type = "gauss",
                                          double size = NA_REAL)
{
  return dust::normalize_1d(std::move(y), type, size);
}
