#pragma once
#include "irred/quantities.hpp"
#include "irred/thermal_neutrino.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <optional>
#include <string_view>

namespace irred::cosmology::detail {
inline constexpr long double thermal_conformal_c_km_s =
    static_cast<long double>(irred::speed_of_light_m_per_s) / 1000;
// Actual retained scalar metadata, captured once by the real consumer after its
// one background preparation. It is not a coefficient/derivative certificate.
struct ThermalRetainedCoefficientWitness {
  numerics::Status status = numerics::Status::invalid_input;
  long double critical_density_ev4 = 0, omega_species_today = 0;
  long double lambda_retained = 0, species_normalization_error = 0;
  double lambda_emitted = 0;
  long double lambda_getter_signed_loss = 0;
  ThermalMomentumMethod momentum_method = ThermalMomentumMethod::direct_adaptive;
};
struct ThermalRetainedCoefficientAccess {
  // Defined in class: implicitly inline. No source vector/P/map/reprepare.
  // Real consumers own once-only capture, matching background/witness copies,
  // move invalidation, actual source/diagnostic work and simultaneous payload.
  static std::optional<ThermalRetainedCoefficientWitness>
  capture(const ThermalBackground &background) noexcept {
    // Prior preparation success does not establish the current rounding mode.
    // Check the retained thermal profile before the one getter conversion.
    if (std::numeric_limits<long double>::digits < 64 ||
        std::numeric_limits<long double>::max_exponent < 16384 ||
        std::fegetround() != FE_TONEAREST ||
        background.status_ != numerics::Status::ok)
      return {};
    const auto emitted = background.omega_lambda();
    if (!emitted || !std::isfinite(*emitted) ||
        !(background.critical_ev4_ > 0) ||
        !std::isfinite(background.critical_ev4_) ||
        !(background.omega_species_ >= 0) ||
        !std::isfinite(background.omega_species_) ||
        !(background.lambda_ >= 0) || !std::isfinite(background.lambda_) ||
        !(background.normalization_error_ >= 0) ||
        !std::isfinite(background.normalization_error_))
      return {};
    ThermalRetainedCoefficientWitness out;
    out.status = numerics::Status::ok;
    out.critical_density_ev4 = background.critical_ev4_;
    out.omega_species_today = background.omega_species_;
    out.lambda_retained = background.lambda_;
    out.species_normalization_error = background.normalization_error_;
    out.lambda_emitted = *emitted;
    out.lambda_getter_signed_loss =
        background.lambda_ - static_cast<long double>(*emitted);
    out.momentum_method = background.method_;
    if (!std::isfinite(out.lambda_getter_signed_loss))
      return {};
    return out;
  }
};
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
  long double p_shadow = 0, p_shadow_n = 0, p_shadow_radius = 0;
  long double actual_p_minus_shadow = 0;
  long double actual_p_shadow_arithmetic_estimate = 0;
  long double actual_p_source_radius = 0;
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
  // For the admitted same-owner source/query chain, P is positive and the
  // supplied cold/radiation/Lambda numerators are nonnegative. Guard every
  // stored derived coordinate and literal diagnostic before publishing ok;
  // algebraic availability does not earn either optional error bundle.
  if (!(out.p > 0) || !(out.hcal > 0) || !std::isfinite(out.p) ||
      !std::isfinite(out.hcal) || !std::isfinite(out.g) ||
      !std::isfinite(out.closure_defect) ||
      !std::isfinite(out.background_identity_defect) ||
      !std::isfinite(out.acceleration_defect))
    out.status = numerics::Status::overflow;
  for (const auto value : {out.p_error, out.fb, out.fc, out.fg, out.fr, out.fl,
                           out.x2, out.cold_fraction, out.enthalpy_fraction,
                           out.lambda_cast_error})
    if (!(value >= 0) || !std::isfinite(value))
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
// Opt-in dependency diagnostic for the physical mapped photon/cold matter
// source only. The same consumer owns the once-captured retained coefficient
// witness, original four mapper records and this already charged epoch/query.
// No species, other radiation, remap, getter, P/H query or alternate H here.
// These conditional estimates inherit the mapper's deterministic operation
// estimates. They are not uncertainty in inputs or a universal certificate.
// The old fluid/FD conversion above does not invoke this function.
inline numerics::Status thermal_conformal_no_species_diagnostics(
    const ThermalBackground &background,
    const ThermalRetainedCoefficientWitness &retained,
    const std::array<ThermalScalarMapWitness, 4> &mapping,
    long double a, ThermalConformalEpoch &epoch) {
  using S = numerics::Status;
  using W = long double;
  epoch.forward.reset(); epoch.shadow.reset();
  const auto &source = background.source();
  if (std::numeric_limits<W>::digits < 64 ||
      std::numeric_limits<W>::max_exponent < 16384 ||
      std::fegetround() != FE_TONEAREST || background.status() != S::ok ||
      retained.status != S::ok || epoch.status != S::ok ||
      !epoch.raw_scaled_query || epoch.raw_scaled_query->status != S::ok ||
      epoch.raw_scaled_query->a4_e2 != epoch.p ||
      epoch.raw_scaled_query->error_estimate != epoch.p_error ||
      retained.momentum_method != background.momentum_method() ||
      !(a > 0 && a <= 1) || !std::isfinite(a))
    return S::conditioning_budget_exceeded;
  if (!source.species.empty() || source.omega_massless_nonphoton != 0 ||
      !(source.omega_gamma > 0) || !(source.omega_b > 0) || !(source.omega_cdm > 0) ||
      retained.omega_species_today != 0 || retained.species_normalization_error != 0 ||
      epoch.fr != epoch.fg || !(epoch.x2 >= std::numeric_limits<W>::min()))
    return S::outside_domain;
  const auto normal = [](W v) {
    return std::isfinite(v) && (v == 0 || std::abs(v) >= std::numeric_limits<W>::min());
  };
  const std::array<W, 4> emitted{source.omega_gamma, source.omega_b,
                                source.omega_cdm, source.omega_massless_nonphoton};
  std::array<W, 4> radius{};
  for (unsigned i = 0; i < 4; ++i) {
    const auto &m = mapping[i];
    if (W(m.emitted_value) != emitted[i] || !(m.wide_value >= 0) ||
        !(m.wide_operation_estimate >= 0) || !(m.measured_absolute_cast_loss >= 0) ||
        !normal(m.wide_value) || !normal(m.wide_operation_estimate) ||
        !normal(m.measured_absolute_cast_loss) ||
        m.measured_absolute_cast_loss != std::abs(m.wide_value - W(m.emitted_value)))
      return S::conditioning_budget_exceeded;
    radius[i] = m.wide_operation_estimate + m.measured_absolute_cast_loss;
    if (!normal(radius[i])) return S::conditioning_budget_exceeded;
  }
  if (!(retained.lambda_retained >= 0) || !normal(retained.lambda_retained) ||
      !(retained.lambda_emitted >= 0) || !normal(retained.lambda_emitted) ||
      !normal(retained.lambda_getter_signed_loss) ||
      retained.lambda_getter_signed_loss != retained.lambda_retained - W(retained.lambda_emitted))
    return S::conditioning_budget_exceeded;
  constexpr W e = std::numeric_limits<W>::epsilon();
  const W a2 = a * a, a4 = a2 * a2;
  if (!(a4 > 0) || !normal(a4)) return S::outside_domain;
  // Largest-first preparation has twenty retained subtractions. Mapping
  // effects enter its closure once; the allowance also covers radius assembly.
  const W lambda_radius = radius[0] + radius[1] + radius[2] + radius[3] +
      32 * e * (1 + emitted[0] + emitted[1] + emitted[2] + emitted[3]);
  const W matter = (emitted[1] + emitted[2]) * a;
  const W vacuum = retained.lambda_retained * a4;
  ThermalConformalShadowDerivativeDiagnostic shadow;
  shadow.p_shadow = emitted[0] + matter + vacuum;
  shadow.p_shadow_n = matter + 4 * vacuum;
  const W p_arithmetic = 64 * e * (emitted[0] + matter + vacuum);
  shadow.p_shadow_radius = radius[0] + a * (radius[1] + radius[2]) +
                           a4 * lambda_radius + p_arithmetic;
  const W pn_radius = a * (radius[1] + radius[2]) + 4 * a4 * lambda_radius +
                     64 * e * (matter + 4 * vacuum);
  shadow.actual_p_minus_shadow = epoch.p - shadow.p_shadow;
  shadow.actual_p_shadow_arithmetic_estimate =
      std::abs(shadow.actual_p_minus_shadow) + p_arithmetic;
  shadow.actual_p_source_radius = std::abs(shadow.actual_p_minus_shadow) +
                                  shadow.p_shadow_radius + epoch.p_error;
  const W shadow_lower = shadow.p_shadow - shadow.p_shadow_radius;
  const W actual_lower = epoch.p - shadow.actual_p_source_radius;
  if (!(shadow_lower > 0) || !(actual_lower > 0) ||
      !normal(shadow_lower) || !normal(actual_lower))
    return S::conditioning_budget_exceeded;
  shadow.q = shadow.p_shadow_n / shadow.p_shadow;
  const W q_product = shadow.q * shadow.p_shadow;
  shadow.q_radius = (pn_radius + std::abs(shadow.q) * shadow.p_shadow_radius +
      std::abs(shadow.p_shadow_n - q_product) +
      8 * e * (std::abs(shadow.p_shadow_n) + std::abs(q_product))) /
      shadow_lower + 8 * e * std::abs(shadow.q);
  shadow.ell = shadow.q / 2 - 1;
  shadow.ell_radius = shadow.q_radius / 2 + 8 * e * (1 + std::abs(shadow.ell));
  ThermalConformalForwardDiagnostic forward;
  const auto fraction_radius = [&](W value, W numerator, W numerator_radius) {
    const W product = value * epoch.p;
    return (numerator_radius + std::abs(value) * shadow.actual_p_source_radius +
            std::abs(numerator - product) +
            8 * e * (std::abs(numerator) + std::abs(product))) /
            actual_lower + 8 * e * std::abs(value);
  };
  forward.fb_absolute_estimate = fraction_radius(epoch.fb, a * emitted[1], a * radius[1]);
  forward.fc_absolute_estimate = fraction_radius(epoch.fc, a * emitted[2], a * radius[2]);
  forward.fg_absolute_estimate = fraction_radius(epoch.fg, emitted[0], radius[0]);
  forward.fr_absolute_estimate = forward.fg_absolute_estimate;
  forward.fl_absolute_estimate = fraction_radius(epoch.fl, W(retained.lambda_emitted) * a4,
      a4 * (lambda_radius + std::abs(retained.lambda_getter_signed_loss)) +
      8 * e * std::abs(W(retained.lambda_emitted) * a4));
  forward.g_absolute_estimate =
      (forward.fb_absolute_estimate + forward.fc_absolute_estimate) / 2 +
      2 * forward.fl_absolute_estimate +
      16 * e * (1 + std::abs(epoch.cold_fraction) / 2 + 2 * std::abs(epoch.fl));
  // The consumer needs the actual used g versus the literal shadow derivative,
  // including their measured center difference and mapped-source variation.
  forward.g_absolute_estimate = std::max(forward.g_absolute_estimate,
      std::abs(epoch.g - shadow.ell) + shadow.ell_radius +
      16 * e * (1 + std::abs(epoch.cold_fraction) / 2 + 2 * std::abs(epoch.fl)));
  forward.closure_assembly_absolute_estimate =
      16 * e * (1 + epoch.cold_fraction + epoch.fr + epoch.fl);
  forward.acceleration_assembly_absolute_estimate =
      16 * e * (std::abs(epoch.g) + 1 + 1.5L * epoch.enthalpy_fraction);
  // Postcheck the ACTUAL H against its squared defining relation. For positive
  // hhat and Hsource, |hhat-Hsource| <= |hhat²-Hsource²|/hhat. Incorporate the
  // source-P envelope inside that numerator: no second sqrt/H and no omitted
  // numerical/source cross-term. Fixed conversion operation effects are named.
  const W conversion = W(source.h0_km_s_mpc) / thermal_conformal_c_km_s / a;
  const W coefficient = conversion * conversion;
  const W actual_h2 = epoch.hcal * epoch.hcal;
  const W target_h2 = coefficient * epoch.p;
  const W source_h2_radius = coefficient * shadow.actual_p_source_radius;
  forward.hcal_absolute_estimate =
      (std::abs(actual_h2 - target_h2) + source_h2_radius +
       256 * e * (std::abs(actual_h2) + std::abs(target_h2) + source_h2_radius)) /
      epoch.hcal;
  if (!(epoch.hcal > forward.hcal_absolute_estimate) ||
      !normal(coefficient) || !normal(actual_h2) || !normal(target_h2) ||
      !normal(source_h2_radius))
    return S::conditioning_budget_exceeded;
  const W relative_h = forward.hcal_absolute_estimate /
                       (epoch.hcal - forward.hcal_absolute_estimate);
  forward.x2_absolute_estimate = std::abs(epoch.x2) *
      (2 * relative_h + relative_h * relative_h) +
      16 * e * std::abs(epoch.x2) * (1 + relative_h) * (1 + relative_h);
  shadow.dp = shadow.q - 2 * (epoch.g + 1);
  shadow.dp_radius = shadow.q_radius + 2 * forward.g_absolute_estimate +
                    8 * e * (std::abs(shadow.q) + 2 * std::abs(epoch.g + 1));
  shadow.source_law_id = "mapped-photon-cold-matter-retained-polynomial/conditional-diagnostic/v1";
  for (W v : {lambda_radius, p_arithmetic, pn_radius, shadow.p_shadow,
              shadow.p_shadow_n, shadow.p_shadow_radius,
              shadow.actual_p_shadow_arithmetic_estimate, shadow.actual_p_source_radius,
              shadow.q_radius, shadow.ell_radius, shadow.dp_radius,
              forward.hcal_absolute_estimate, forward.x2_absolute_estimate,
              forward.fb_absolute_estimate, forward.fc_absolute_estimate,
              forward.fg_absolute_estimate, forward.fr_absolute_estimate,
              forward.fl_absolute_estimate, forward.g_absolute_estimate,
              forward.closure_assembly_absolute_estimate,
              forward.acceleration_assembly_absolute_estimate})
    if (!(v >= 0) || !normal(v)) return S::conditioning_budget_exceeded;
  for (W v : {shadow.q, shadow.ell, shadow.dp, shadow.actual_p_minus_shadow})
    if (!normal(v)) return S::conditioning_budget_exceeded;
  epoch.forward = forward; epoch.shadow = shadow;
  return S::ok;
}
} // namespace irred::cosmology::detail
