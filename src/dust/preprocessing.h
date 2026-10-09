/// Data preprocessing: noise level and normalization for the default penalty

#ifndef DUST_PREPROCESSING_H
#define DUST_PREPROCESSING_H

#include <string>
#include <vector>

namespace dust {

/// standard deviation of the noise from the differences of y (robust to
/// changes in mean): method "HALL", "MAD" or "SD"
double sd_diff(const std::vector<double>& y, const std::string& method);

/// y normalized for the model `type` (size: number of trials or successes
/// for "binom" and "negbin")
std::vector<double> normalize_1d(std::vector<double> y, const std::string& type, double size);

} // namespace dust

#endif
