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
  double time_seconds = 0, distance_mpc = 0;
};
inline FlatScale prepare_flat_scale(double value) {
  FlatScale out;
  const auto h = convert({value, Unit::km_per_s_per_mpc, Role::expansion_rate},
                         {Unit::inverse_second, Role::expansion_rate});
  if (h.status != QuantityStatus::ok) {
    out.status = Status::numerical_failure;
    return out;
  }
  const auto time = 1.L / h.target.value;
  const auto length = (long double)speed_of_light_m_per_s * time;
  if (!representable(time) || !representable(length)) {
    out.status = Status::numerical_failure;
    return out;
  }
  const auto d = convert({(double)length, Unit::metre, Role::physical_length,
                          Frame::none, LengthConvention::physical},
                         {Unit::megaparsec, Role::physical_length, Frame::none,
                          LengthConvention::physical});
  if (d.status != QuantityStatus::ok) {
    out.status = Status::numerical_failure;
    return out;
  }
  out.time_seconds = (double)time;
  out.distance_mpc = d.target.value;
  out.status = Status::ok;
  return out;
}
} // namespace irred::cosmology::detail
