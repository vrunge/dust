#ifndef DUST_1D_VALIDATION_H
#define DUST_1D_VALIDATION_H

#include <Rcpp.h>
#include <cmath>

namespace dust_validation {

inline double penalty(Rcpp::Nullable<double> supplied, double current,
                      std::size_t first_batch_size, bool initialized)
{
  if (supplied.isNull())
    return initialized ? current : 2.0 * std::log(static_cast<double>(first_batch_size));
  const double value = Rcpp::as<double>(supplied);
  if (!std::isfinite(value) || value < 0.0)
    Rcpp::stop("penalty must be a finite nonnegative number");
  if (initialized && value != current)
    Rcpp::stop("penalty cannot change after the first nonempty append");
  return value;
}

} // namespace dust_validation

#endif
