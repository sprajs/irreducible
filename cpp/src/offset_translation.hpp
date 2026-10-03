#pragma once
#include "irred/numerics.hpp"
#include <cmath>
#include <limits>
#include <span>
#include <vector>
namespace irred::detail {
// Exact translation is a declared domain, not an omitted arithmetic error.
inline numerics::Status center_exact(std::span<const double> y,
                                     std::span<const double> offsets,
                                     std::vector<double> &out) {
  if (y.size() != offsets.size())
    return numerics::Status::invalid_input;
  out.resize(y.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    if (!std::isfinite(y[i]) || (y[i] != 0 && !std::isnormal(y[i])))
      return numerics::Status::nonfinite_input;
    const long double exact = static_cast<long double>(y[i]) - offsets[i];
    const double rounded = static_cast<double>(exact);
    if (!std::isfinite(rounded) || (rounded != 0 && !std::isnormal(rounded)) ||
        (exact != 0 && rounded == 0))
      return numerics::Status::outside_domain;
    const double neg = -offsets[i], z = rounded - y[i];
    const double low = (y[i] - (rounded - z)) + (neg - z);
    if (low != 0 || static_cast<long double>(rounded) != exact)
      return numerics::Status::conditioning_budget_exceeded;
    out[i] = rounded;
  }
  return numerics::Status::ok;
}
// Translate a derived conditional mean, propagating actual storage rounding.
inline numerics::Status
add_offsets(std::span<const double> mean, std::span<const double> errors,
            std::span<const double> offsets, double sensitivity,
            std::vector<double> &values, std::vector<double> &reported_errors) {
  if (mean.size() != errors.size() || mean.size() != offsets.size())
    return numerics::Status::invalid_input;
  values.resize(mean.size());
  reported_errors.resize(mean.size());
  for (std::size_t i = 0; i < mean.size(); ++i) {
    const long double wide = static_cast<long double>(mean[i]) + offsets[i];
    const double value = static_cast<double>(wide);
    const long double err = errors[i] + std::abs(wide - value) +
                            4 * std::numeric_limits<long double>::epsilon() *
                                (std::abs(mean[i]) + std::abs(offsets[i]));
    double error = static_cast<double>(err);
    if (static_cast<long double>(error) < err)
      error = std::nextafter(error, std::numeric_limits<double>::infinity());
    if (!std::isfinite(value) || (value != 0 && !std::isnormal(value)) ||
        (wide != 0 && value == 0) || !std::isfinite(error) ||
        (error != 0 && !std::isnormal(error)) || (err > 0 && error == 0))
      return numerics::Status::outside_domain;
    if (err > sensitivity * (1 + std::abs(wide)))
      return numerics::Status::conditioning_budget_exceeded;
    values[i] = value;
    reported_errors[i] = error;
  }
  return numerics::Status::ok;
}
} // namespace irred::detail
