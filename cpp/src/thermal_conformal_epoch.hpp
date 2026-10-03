#pragma once
#include "irred/quantities.hpp"
#include "irred/thermal_neutrino.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string_view>

namespace irred::cosmology::detail {
inline constexpr long double thermal_conformal_c_km_s =
    static_cast<long double>(irred::speed_of_light_m_per_s) / 1000;
// These complete diagnostics are unavailable until their source/arithmetic
// ownership is earned. Entries are absolute estimates/radii, in the same units
// as their epoch coordinates (hcal is Mpc^-1; the rest are dimensionless).
// Every entry must be finite and nonnegative; positive entries remain normal.
struct ThermalConformalForwardDiagnostic {
  long double hcal_absolute_estimate = 0, x2_absolute_estimate = 0;
  long double fb_absolute_estimate = 0, fc_absolute_estimate = 0;
  long double fg_absolute_estimate = 0, fr_absolute_estimate = 0;
  long double fl_absolute_estimate = 0, g_absolute_estimate = 0;
  // Supplementary loss in the literal diagnostic assembly, not another force.
  long double closure_assembly_absolute_estimate = 0;
  long double acceleration_assembly_absolute_estimate = 0;
};
struct ThermalConformalShadowDerivativeDiagnostic {
  // Future source semantics require q=P_shadow_N/P_shadow and an earned
  // P_actual/P_shadow discrepancy envelope. These are never derivatives of a
  // rounded callback or a quotient with P_actual as the shadow denominator.
  long double q = 0, ell = 0, dp = 0;
  long double q_radius = 0, ell_radius = 0, dp_radius = 0;
  std::string_view source_law_id; // Must have static lifetime.
};
// Shared conversion of the actual retained thermal P(a) to conformal length
// units. radiation is the caller's explicitly identified total radiation
// numerator: supplied perfect fluid, explicit massless FD and ideal photons
// remain distinct sources. fg is its photon subrole and is never added again.
struct ThermalConformalEpoch {
  numerics::Status status = numerics::Status::invalid_input;
  // status==ok is the availability discriminator for these derived values.
  // The independent raw query receipt survives either query/conversion refusal.
  long double p = 0, p_error = 0, hcal = 0;
  long double fb = 0, fc = 0, fg = 0, fr = 0, fl = 0, x2 = 0, g = 0;
  long double cold_fraction = 0, enthalpy_fraction = 0;
  long double closure_defect = 0, background_identity_defect = 0;
  long double acceleration_defect = 0, lambda_cast_error = 0;
  std::size_t momentum_callbacks = 0;
  std::optional<ThermalScaledExpansion> raw_scaled_query;
  // Neither complete bundle is earned by this algebraic conversion. Consumers
  // requiring them must refuse their own dependency admission if absent.
  std::optional<ThermalConformalForwardDiagnostic> forward;
  std::optional<ThermalConformalShadowDerivativeDiagnostic> shadow;
};
// Private borrowed import: scaled must be the single already charged query of
// this same background, a and retained momentum method, under the same frozen
// arithmetic profile. ThermalScaledExpansion carries no owner/a/method tags;
// this is an internal call-chain contract, never verified by a second query.
// Import performs no query or callback work. Its callback receipt is not a
// second charge. Consumers keep dependency refusal separate from epoch status.
inline ThermalConformalEpoch thermal_conformal_epoch_from_scaled(
    const ThermalBackground &background, const ThermalScaledExpansion &scaled,
    long double a, long double k, long double radiation) {
  ThermalConformalEpoch out;
  out.raw_scaled_query = scaled;
  out.momentum_callbacks = scaled.callbacks;
  if (!(a > 0) || !(k >= 0) || !(radiation >= 0) ||
      !std::isfinite(a) || !std::isfinite(k) || !std::isfinite(radiation))
    return out;
  out.status = scaled.status;
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
  out.fb = static_cast<long double>(source.omega_b) * a / out.p;
  out.fc = static_cast<long double>(source.omega_cdm) * a / out.p;
  out.fg = static_cast<long double>(source.omega_gamma) / out.p;
  out.fr = radiation / out.p;
  out.fl = static_cast<long double>(*lambda) * a4 / out.p;
  out.x2 = (k / out.hcal) * (k / out.hcal);
  out.cold_fraction = out.fb + out.fc;
  out.enthalpy_fraction = out.cold_fraction + 4 * out.fr / 3;
  out.g = -1 + out.cold_fraction / 2 + 2 * out.fl;
  // Preserve the original grouping at fb==0. The separately rounded consumer
  // Rbg assembly below need not be identical to this legacy identity diagnostic.
  out.closure_defect = out.cold_fraction + out.fr + out.fl - 1;
  out.background_identity_defect =
      out.g - 1 + 1.5L * out.cold_fraction + 2 * out.fr;
  out.acceleration_defect = out.g - 1 + 1.5L * out.enthalpy_fraction;
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
inline ThermalConformalEpoch thermal_conformal_epoch(
    const ThermalBackground &background, long double a, long double k,
    long double radiation, ThermalPolicy policy) {
  // An invalid wrapper preflight performs no query and has no raw receipt.
  if (!(a > 0) || !(k >= 0) || !(radiation >= 0) ||
      !std::isfinite(a) || !std::isfinite(k) || !std::isfinite(radiation))
    return {};
  const auto scaled = background.scaled_expansion(a, policy);
  return thermal_conformal_epoch_from_scaled(background, scaled, a, k, radiation);
}
} // namespace irred::cosmology::detail
