#pragma once
// Shared flat geometry/quantity arithmetic. Extraction preserves legacy operand
// and cast order; providers supply expansion/integrals and derivative
// semantics.
#include "irred/background.hpp"
#include <cmath>
#include <limits>
namespace irred::cosmology::detail {
inline bool physical_representable(long double v) {
  return std::isfinite(v) &&
         std::abs(v) <= std::numeric_limits<double>::max() &&
         (v == 0 || std::abs(v) >= std::numeric_limits<double>::min());
}
inline bool representable(long double v) {
  return std::isfinite(v) && std::abs(v) <= std::numeric_limits<double>::max();
}
struct FlatScale {
  Status status = Status::invalid_input;
  double time_seconds = 0, distance_mpc = 0;
};
inline FlatScale prepare_flat_scale(double h0_value) {
  FlatScale out;
  const auto h0 =
      convert({h0_value, Unit::km_per_s_per_mpc, Role::expansion_rate},
              {Unit::inverse_second, Role::expansion_rate});
  if (h0.status != QuantityStatus::ok) {
    out.status = Status::numerical_failure;
    return out;
  }
  const auto time = 1.L / h0.target.value;
  const auto length = static_cast<long double>(speed_of_light_m_per_s) * time;
  if (!representable(time) || !representable(length)) {
    out.status = Status::numerical_failure;
    return out;
  }
  const auto dh =
      convert({static_cast<double>(length), Unit::metre, Role::physical_length,
               Frame::none, LengthConvention::physical},
              {Unit::megaparsec, Role::physical_length, Frame::none,
               LengthConvention::physical});
  if (dh.status != QuantityStatus::ok) {
    out.status = Status::numerical_failure;
    return out;
  }
  out.time_seconds = static_cast<double>(time);
  out.distance_mpc = dh.target.value;
  out.status = Status::ok;
  return out;
}
struct FlatGeometry {
  Status status = Status::numerical_failure;
  numerics::Status numerical_status = numerics::Status::outside_domain;
  long double dc = 0, da = 0, shape = 0, dl = 0, lookback = 0, volume = 0;
};
inline FlatGeometry flat_geometry(const Query &query, long double e,
                                  double radial, double clock,
                                  double distance_mpc, double time_seconds) {
  FlatGeometry out;
  const auto u = 1 + static_cast<long double>(query.z_expansion);
  const auto dc = static_cast<long double>(distance_mpc) * radial;
  const auto da = dc / u;
  const auto shape = (1 + static_cast<long double>(query.z_observer)) * radial;
  const auto dl = static_cast<long double>(distance_mpc) * shape;
  const auto lookback = static_cast<long double>(time_seconds) * clock;
  const auto volume =
      static_cast<long double>(distance_mpc) * dc * dc / e; // dVc/(dz dOmega)
  if (query.z_expansion > 0 &&
      (!(dc > 0) || !(da > 0) || !(shape > 0) || !(dl > 0) || !(lookback > 0) ||
       !(volume > 0))) {
    return out;
  }
  if (!physical_representable(dc) || !physical_representable(da) ||
      !representable(shape) || !physical_representable(dl) ||
      !physical_representable(lookback) || !physical_representable(volume)) {
    for (auto x : {dc, da, shape, dl, lookback, volume})
      if (!std::isfinite(x) || std::abs(x) > std::numeric_limits<double>::max())
        out.numerical_status = numerics::Status::overflow;
    return out;
  }
  out.status = Status::ok;
  out.dc = dc;
  out.da = da;
  out.shape = shape;
  out.dl = dl;
  out.lookback = lookback;
  out.volume = volume;
  return out;
}
} // namespace irred::cosmology::detail
