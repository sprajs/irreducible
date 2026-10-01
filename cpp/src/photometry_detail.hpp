#pragma once
#include "irred/photometry.hpp"
#include "irred/quantities.hpp"
#include <cmath>
#include <limits>
#include <numbers>
namespace irred::photometry::detail {
using Wide = long double;
struct Interval { Wide lower, upper; };
// Error-free TwoSum under the declared round-to-nearest contract.
inline Interval add(Wide a, Wide b) noexcept {
  const Wide value = a + b;
  const Wide bb = value - a;
  const Wide error = (a - (value - bb)) + (b - bb);
  return {error < 0 ? std::nextafter(value, -std::numeric_limits<Wide>::infinity()) : value,
          error > 0 ? std::nextafter(value, std::numeric_limits<Wide>::infinity()) : value};
}
inline Interval multiply(Wide a, Wide b) noexcept {
  const Wide value = a * b;
  const Wide error = std::fma(a, b, -value);
  return {error < 0 ? std::nextafter(value, -std::numeric_limits<Wide>::infinity()) : value,
          error > 0 ? std::nextafter(value, std::numeric_limits<Wide>::infinity()) : value};
}
inline Interval scale(Interval a, Wide b) noexcept {
  return {multiply(a.lower, b).lower, multiply(a.upper, b).upper};
}
inline void fail(Outcome &out, numerics::Status cause) noexcept {
  out.availability = Availability::failed;
  out.numerical_status = cause;
  out.value.reset();
}
inline void store(Outcome &out, Wide value) noexcept {
  if (!std::isfinite(value) || value > std::numeric_limits<double>::max()) {
    fail(out, numerics::Status::overflow); return;
  }
  const double rounded = static_cast<double>(value);
  if (value > 0 && rounded == 0) {
    fail(out, numerics::Status::outside_domain); return;
  }
  if (value > 0 && std::abs(static_cast<Wide>(rounded) - value) / value > 1e-12L) {
    fail(out, numerics::Status::conditioning_budget_exceeded); return;
  }
  out.availability = Availability::available;
  out.numerical_status = numerics::Status::ok;
  out.value = rounded;
}
inline Wide spectral_flux(Wide luminosity, Wide distance, Wide redshift_factor) noexcept {
  return luminosity / (4 * std::numbers::pi_v<Wide> * distance * distance * redshift_factor);
}
} // namespace irred::photometry::detail
