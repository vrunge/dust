#ifndef DUST_PREPROCESSING_H
#define DUST_PREPROCESSING_H

#include <string>
#include <vector>

namespace dust {

/// noise level from the differences of y: "HALL", "MAD" or "SD"
double sd_diff(const std::vector<double>& y, const std::string& method);

/// size: number of trials (binom) or successes (negbin)
std::vector<double> normalize_1d(std::vector<double> y, const std::string& type, double size);

} // namespace dust

#endif
