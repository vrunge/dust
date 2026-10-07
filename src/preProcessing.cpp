#include <Rcpp.h>

using namespace Rcpp;

#include <Rcpp.h>
#include <cmath>
#include <algorithm> // for std::max_element, for std::all_of and std::isfinite
#include <numeric>

#include "preProcessing.h"


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

//' data_normalization_1D
//'
//' @name data_normalization_1D
//'
//' @description Normalization of the data before using dust.1D with the default penalty
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
