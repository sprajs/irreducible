#pragma once
#include "irred/background.hpp"
#include <cmath>
#include <limits>
namespace irred::cosmology::detail {
inline bool representable(long double v) {
  return std::isfinite(v) && std::abs(v) <= std::numeric_limits<double>::max();
}
inline bool physical_representable(long double v) {
  return representable(v) &&
         (v == 0 || std::abs(v) >= std::numeric_limits<double>::min());
}
struct FlatScale {
  Status status = Status::invalid_input;
  long double time_seconds = 0, distance_mpc = 0;
};
inline FlatScale prepare_flat_scale(double value) {
  FlatScale out;
  if (!std::isfinite(value) || !(value > 0))
    return out;
  // H0 is supplied in km/s/Mpc. Keep the unit scales wide and validate only
  // requested final seconds/Mpc/volume values in the consuming projection.
  out.time_seconds = (megaparsec_in_metres_wide() / 1000) / value;
  out.distance_mpc = ((long double)speed_of_light_m_per_s / 1000) / value;
  out.status = Status::ok;
  return out;
}
} // namespace irred::cosmology::detail
