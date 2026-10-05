#include <Rcpp.h>

using namespace Rcpp;

#include <Rcpp.h>
#include <cmath>
#include <algorithm> // for std::max_element, for std::all_of and std::isfinite
#include <numeric>

#include "preProcessing.h"


//' Calculate Standard Deviation or MAD of Differences in a Numeric Vector
//'
//' The `sdDiff` function calculates a measure of variability (standard deviation or MAD)
//' of a numeric vector.
//' It supports three methods: "HALL", "MAD", and "SD".
//'
//' @param y A numeric vector.
//' @param method A character string specifying the method to use.
//'   Options are: \code{"HALL"}, \code{"MAD"}, and \code{"SD"}.
//'   Default is \code{"HALL"}.
//'
//' @return A numeric value representing the calculated measure of variability
//'   according to the specified method.
//'   \itemize{
//'     \item \code{"HALL"}: Calculates the standard deviation using a specific weighted
//'           difference method (HALL method).
//'     \item \code{"MAD"}: Returns the MAD (Median Absolute Deviation) of the differences
//'           between consecutive elements.
//'     \item \code{"SD"}: Returns the standard deviation of the differences
//'           between consecutive elements.
//'   }
//'
//' @examples
//' ### 3 segments of 300 points, mean shifts but constant noise sd = 2:
//' ### all three estimators should recover a value close to 2.
//' set.seed(30)
//' y <- dataGenerator_1D(chpts = c(300, 600, 900), parameters = c(0, 1, 0),
//'                        sdNoise = 2, type = "gauss")
//' sdDiff(y, "HALL")
//' sdDiff(y, "MAD")
//' sdDiff(y, "SD")
//'
//' @export
// [[Rcpp::export]]
double sdDiff(std::vector<double>& y, std::string method = "HALL")
{
  ///////////////////////  HALL
  ///////////////////////  HALL
  ///////////////////////  HALL
  if(method == "HALL")
  {
    if(!std::all_of(y.begin(), y.end(), [](double val) { return std::isfinite(val); }) || y.size() < 5)
    {
      Rcpp::stop("y is not a numeric vector or length < 5 (the HALL method cannot be used)");
    }

    int n = y.size();
    Rcpp::NumericVector wei = {0.1942, 0.2809, 0.3832, -0.8582};
    Rcpp::NumericMatrix mat(4, n);

    // Constructing the matrix `mat`
    for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < n; j++)
      {
        mat(i, j) = wei[i] * y[j];
      }
    }

    // Adjusting the elements according to the R function
    for(int j = 0; j < n - 1; j++){mat(1, j) = mat(1, j + 1);}
    for(int j = 0; j < n - 2; j++){mat(2, j) = mat(2, j + 2);}
    for(int j = 0; j < n - 3; j++){mat(3, j) = mat(3, j + 3);}

    // Computing the result
    double sumSquares = 0.0;
    double columnSum = 0.0;
    for(int j = 0; j < n - 3; j++)
    {
      columnSum = 0.0;
      for(int i = 0; i < 4; i++){columnSum += mat(i, j);}
      sumSquares += columnSum * columnSum;
    }
    return std::sqrt(sumSquares / (n - 3));
  }
  ///////////////////////  MAD
  ///////////////////////  MAD
  ///////////////////////  MAD
  if(method == "MAD")
  {
    if(!std::all_of(y.begin(), y.end(), [](double val) { return std::isfinite(val); }) || y.size() < 2)
    {
      Rcpp::stop("y is not a numeric vector or length < 2 (the MAD method cannot be used)");
    }
    int n = y.size() - 1;

    Rcpp::NumericVector result(n);
    double sqrt2 = std::sqrt(2.0);
    for (int i = 0; i < n ; ++i) {result[i] = (y[i + 1] - y[i]) / sqrt2;}

    // Calculate the median
    std::vector<double> sorted_x(result.begin(), result.end());
    std::nth_element(sorted_x.begin(), sorted_x.begin() + n / 2, sorted_x.end());
    double median_x = sorted_x[n / 2];
    if (n % 2 == 0)
    {
      std::nth_element(sorted_x.begin(), sorted_x.begin() + n / 2 - 1, sorted_x.end());
      median_x = (median_x + sorted_x[n / 2 - 1]) / 2.0;
    }

    // Calculate the absolute deviations from the median
    std::vector<double> abs_dev(n);
    for (int i = 0; i < n; i++)
    {
      abs_dev[i] = std::abs(result[i] - median_x);
    }

    // Calculate the median of the absolute deviations
    std::nth_element(abs_dev.begin(), abs_dev.begin() + n / 2, abs_dev.end());
    double mad = abs_dev[n / 2];
    if (n % 2 == 0)
    {
      std::nth_element(abs_dev.begin(), abs_dev.begin() + n / 2 - 1, abs_dev.end());
      mad = (mad + abs_dev[n / 2 - 1]) / 2.0;
    }
    return mad * 1.4826;
  }

  ///////////////////////  SD
  ///////////////////////  SD
  ///////////////////////  SD
  if(method == "SD")
  {
    int n = y.size() - 1;
    if(!std::all_of(y.begin(), y.end(), [](double val) { return std::isfinite(val); }) || y.size() < 2)
    {
      Rcpp::stop("y is not a numeric vector or length < 2 (the SD method cannot be used)");
    }

    Rcpp::NumericVector diff_y(n);
    double sqrt2 = std::sqrt(2.0);
    for (int i = 0; i < n; ++i) {diff_y[i] = (y[i + 1] - y[i]) / sqrt2;}

    // Calculate the standard deviation
    double mean = Rcpp::mean(diff_y);
    double sum_sq_diff = 0.0;
    for (int i = 0; i < n; ++i){sum_sq_diff += std::pow(diff_y[i] - mean, 2);}

    double variance = sum_sq_diff / (n - 1);
    return std::sqrt(variance);
  }
  Rcpp::stop("method must be one of HALL, MAD, or SD");
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

//' Data Normalization Function
//'
//' @name data_normalization_1D
//'
//' @description
//' Normalizes the input time series data `y` according to the specified `type`.
//' The normalization process depends on the statistical model type, which can be one of the following:
//' "gauss" (Gaussian/normal distribution), "exp" (exponential distribution),
//' "poisson" (Poisson distribution), "geom" (geometric distribution),
//' "bern" (Bernoulli distribution), "binom" (binomial distribution),
//' "negbin" (negative binomial distribution), or "variance"
//' Poisson scaling changes the likelihood scale; supply a penalty appropriate
//' for the scaled data. Binomial and Negative Binomial normalization use a
//' known size rather than estimating it from the observations.
//'
//' @param y A numeric vector representing the time series to be normalized and then segmented.
//' @param type A string specifying the model type for normalization.
//' The available options are "gauss", "exp", "poisson", "geom", "bern", "binom", "negbin", "variance".
//' The default is "gauss".
//' @param size Known number of Binomial trials or Negative Binomial size.
//'   Required for those two models; the observed maximum and an estimated
//'   dispersion are not substitutes for the known model parameter.
//' @return A numeric vector that is the normalized version of the input time series `y`.
//' @examples
//' ### Gaussian: 3 segments of 300 points, noise sd = 2 -- normalization
//' ### rescales by an sdDiff() estimate of the noise, so the result has
//' ### noise sd close to 1.
//' set.seed(40)
//' y <- dataGenerator_1D(chpts = c(300, 600, 900), parameters = c(0, 1, 0),
//'                        sdNoise = 2, type = "gauss")
//' y_norm <- data_normalization_1D(y, type = "gauss")
//' sd(diff(y_norm)) / sqrt(2)  # close to 1
//'
//' ### Poisson: 3 segments of 300 points, rates 2/8/4 -- normalization
//' ### rescales by the overall mean, so the result has mean close to 1.
//' set.seed(41)
//' y <- dataGenerator_1D(chpts = c(300, 600, 900), parameters = c(2, 8, 4), type = "poisson")
//' y_norm <- data_normalization_1D(y, type = "poisson")
//' mean(y_norm)  # close to 1
//'
//' ### Exponential: 3 segments of 300 points, rates 2/0.5/3
//' set.seed(42)
//' y <- dataGenerator_1D(chpts = c(300, 600, 900), parameters = c(2, 0.5, 3), type = "exp")
//' y_norm <- data_normalization_1D(y, type = "exp")
//' mean(y_norm)  # close to 1
//'
//' @export
// [[Rcpp::export]]
std::vector<double> data_normalization_1D(std::vector<double>& y,
                                          std::string type = "gauss",
                                          double size = NA_REAL)
{
  int n = y.size();
  if (n == 0 || !std::all_of(y.begin(), y.end(), [](double z) { return std::isfinite(z); }))
    Rcpp::stop("y must be a nonempty finite numeric vector");

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if (type == "variance")
  {
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
    for (int i = 0; i < n; ++i){y[i] = y[i] - mean_y;}
    if (std::any_of(y.begin(), y.end(), [](double z) { return z == 0.0; }))
      Rcpp::stop("variance normalization produced a zero residual");
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if (type == "gauss")
  {
    double sdNoise = sdDiff(y);
    if (!(sdNoise > 0.0) || !std::isfinite(sdNoise))
      Rcpp::stop("Gaussian normalization requires positive finite noise scale");
    for (int i = 0; i < n; ++i){y[i] = y[i] / sdNoise;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if(type == "poisson")
  {
    for(int i = 0; i < n; i++){if(y[i] < 0){throw std::range_error("negative data not compatible with poisson model");}}
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
    if (mean_y == 0.0) return y;
    for(int i = 0; i < n; ++i){y[i] = y[i] / mean_y;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if(type == "exp")
  {
    for(int i = 0; i < n; i++){if(y[i] <= 0){throw std::range_error("Exponential data must be strictly positive");}}
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
    for(int i = 0; i < n; ++i) {y[i] = y[i] / mean_y;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if(type == "binom")
  {
    if (!std::isfinite(size) || size <= 0 || std::floor(size) != size)
      Rcpp::stop("Binomial normalization requires a positive integer size (number of trials)");
    for (double value : y)
      if (value < 0 || value > size || std::floor(value) != value)
        Rcpp::stop("Binomial observations must be integer counts between 0 and size");
    for(int i = 0; i < n; ++i) {y[i] = y[i] / size;}
    return y;
  }

  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  //////////  //////////  //////////  //////////
  if(type == "negbin")
  {
    if (!std::isfinite(size) || size <= 0)
      Rcpp::stop("Negative Binomial normalization requires a positive known size");
    for (double value : y)
      if (value < 0 || std::floor(value) != value)
        Rcpp::stop("Negative Binomial observations must be nonnegative integer counts");
    for(size_t i = 0; i < y.size(); i++){y[i] = y[i]/size;}
    return y;
  }

  if(type == "geom")
  {
    if (std::any_of(y.begin(), y.end(), [](double z) { return z < 1.0; }))
      Rcpp::stop("Geometric observations must be at least one");
    return y;
  }

  if(type == "bern")
  {
    if (std::any_of(y.begin(), y.end(), [](double z) { return z < 0.0 || z > 1.0; }))
      Rcpp::stop("Bernoulli observations must be in [0,1]");
    return y;
  }

  Rcpp::stop("Unsupported type specified.");
}
