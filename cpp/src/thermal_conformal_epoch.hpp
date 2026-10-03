#pragma once
#include "irred/quantities.hpp"
#include "irred/thermal_neutrino.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace irred::cosmology::detail {
inline constexpr long double thermal_conformal_c_km_s =
    static_cast<long double>(irred::speed_of_light_m_per_s) / 1000;
// Shared conversion of the retained thermal P(a) to conformal length units.
// The caller supplies its explicitly identified massless radiation fraction:
// supplied perfect fluid and explicit FD radiation are different sources.
struct ThermalConformalEpoch {
  numerics::Status status = numerics::Status::invalid_input;
  long double p = 0, p_error = 0, hcal = 0;
  long double fc = 0, fr = 0, fl = 0, x2 = 0, g = 0;
  long double closure_defect = 0, background_identity_defect = 0;
  long double lambda_cast_error = 0;
  std::size_t momentum_callbacks = 0;
};
inline ThermalConformalEpoch thermal_conformal_epoch(
    const ThermalBackground &background, long double a, long double k,
    long double radiation, ThermalPolicy policy) {
  ThermalConformalEpoch out;
  if (!(a > 0) || !(k >= 0) || !(radiation >= 0) ||
      !std::isfinite(a) || !std::isfinite(k) || !std::isfinite(radiation))
    return out;
  const auto scaled = background.scaled_expansion(a, policy);
  out.status = scaled.status;
  out.momentum_callbacks = scaled.callbacks;
  if (scaled.status != numerics::Status::ok)
    return out;
  const auto lambda = background.omega_lambda();
  if (!lambda) {
    out.status = numerics::Status::invalid_input;
    return out;
  }
  const auto &source = background.source();
  const long double a2 = a * a, a4 = a2 * a2;
  out.p = scaled.a4_e2;
  out.p_error = scaled.error_estimate;
  out.hcal = static_cast<long double>(source.h0_km_s_mpc) /
             thermal_conformal_c_km_s *
             std::sqrt(out.p) / a;
  out.fc = static_cast<long double>(source.omega_cdm) * a / out.p;
  out.fr = radiation / out.p;
  out.fl = static_cast<long double>(*lambda) * a4 / out.p;
  out.x2 = (k / out.hcal) * (k / out.hcal);
  out.g = -1 + out.fc / 2 + 2 * out.fl;
  out.closure_defect = out.fc + out.fr + out.fl - 1;
  out.background_identity_defect = out.g - 1 + 1.5L * out.fc + 2 * out.fr;
  if (*lambda != 0) {
    const double above = std::nextafter(*lambda,
                                      std::numeric_limits<double>::infinity());
    const double below = std::nextafter(*lambda, 0.0);
    out.lambda_cast_error =
        std::max(static_cast<long double>(above) - *lambda,
                 static_cast<long double>(*lambda) - below) / 2;
  }
  if (!(out.hcal > 0) || !std::isfinite(out.hcal) ||
      !std::isfinite(out.x2) || !std::isfinite(out.g))
    out.status = numerics::Status::overflow;
  return out;
}
} // namespace irred::cosmology::detail
