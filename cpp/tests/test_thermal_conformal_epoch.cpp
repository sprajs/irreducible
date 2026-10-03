#include "irred/linear_transfer.hpp"
#include "../src/thermal_conformal_epoch.hpp"
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>

namespace {
using namespace irred::cosmology;
using S = irred::numerics::Status;
using W = long double;
void need(bool value, const char *message) {
  if (!value) {
    std::cerr << "FAIL " << message << '\n';
    std::exit(1);
  }
}
void raw_matches(const detail::ThermalConformalEpoch &epoch,
                 const ThermalScaledExpansion &query) {
  need(epoch.raw_scaled_query &&
           epoch.raw_scaled_query->status == query.status &&
           epoch.raw_scaled_query->a4_e2 == query.a4_e2 &&
           epoch.raw_scaled_query->error_estimate == query.error_estimate &&
           epoch.raw_scaled_query->callbacks == query.callbacks &&
           epoch.momentum_callbacks == query.callbacks,
       "actual raw query and callbacks retained unchanged");
  need(!epoch.forward && !epoch.shadow,
       "unearned complete forward and shadow diagnostics stay absent");
}
void dyadic_source() {
  // Independent exact dyadic control. At a=1/2, Lambda=1/4 and
  // P=7/32+(17/32)/2+(1/4)/16=1/2. No computed expected P/closure is used.
  const ThermalFlatModel source{70, 3. / 32, 1. / 8, 1. / 4, 9. / 32, {}};
  const auto background = prepare_thermal_background(source);
  need(background.status() == S::ok, "dyadic thermal owner");
  ThermalPolicy policy;
  policy.momentum_method = background.momentum_method();
  const auto query = background.scaled_expansion(.5L, policy);
  need(query.status == S::ok && query.a4_e2 == .5L && query.callbacks == 0,
       "actual retained dyadic P");
  const auto epoch = detail::thermal_conformal_epoch_from_scaled(
      background, query, .5L, .001L, 7.L / 32);
  need(epoch.status == S::ok && epoch.p == .5L &&
           epoch.fb == 1.L / 4 && epoch.fc == 9.L / 32 &&
           epoch.fg == 3.L / 16 && epoch.fr == 7.L / 16 &&
           epoch.fl == 1.L / 32 && epoch.cold_fraction == 17.L / 32 &&
           epoch.g == -43.L / 64 && epoch.closure_defect == 0 &&
           epoch.background_identity_defect == 0,
       "independent total cold matter and source-resolved dyadic fractions");
  // Real-algebra Rbg=0. The new grouping includes the rounded 7/12 term;
  // eight wide epsilons cover its finite elementary assembly here. This is a
  // single arithmetic control, not a source/libm or forward-error certificate.
  need(std::abs(epoch.acceleration_defect) <=
           8 * std::numeric_limits<W>::epsilon(),
       "separately assembled acceleration has the exact real zero limit");
  raw_matches(epoch, query);
  const auto reused = detail::thermal_conformal_epoch_from_scaled(
      background, query, .5L, .001L, 7.L / 32);
  raw_matches(reused, query);
  need(reused.status == S::ok && reused.hcal == epoch.hcal &&
           reused.x2 == epoch.x2 && reused.g == epoch.g,
       "paired states import the same one-query stage");
  const auto wrapped = detail::thermal_conformal_epoch(
      background, .5L, .001L, 7.L / 32, policy);
  raw_matches(wrapped, query);
  need(wrapped.status == S::ok && wrapped.hcal == epoch.hcal &&
           wrapped.x2 == epoch.x2 && wrapped.acceleration_defect ==
                                            epoch.acceleration_defect,
       "one-query wrapper and borrowed import share the conversion");
  // Structural invalid-numerator adversary only: retain the REAL same-owner P,
  // but pass a finite unsupported radiation numerator whose quotient overflows.
  // This is no physical comparison and supplies no source/forward certificate.
  const auto overflow = detail::thermal_conformal_epoch_from_scaled(
      background, query, .5L, .001L, std::numeric_limits<W>::max());
  need(overflow.status == S::overflow,
       "nonfinite derived radiation cannot publish algebraic success");
  raw_matches(overflow, query);
  need(prepare_perfect_fluid_transfer(background, 1e-9).status() ==
           S::outside_domain,
       "old public fluid still refuses cold baryons");

  // The same total cold numerator with zero baryons gives the ORIGINAL
  // perfect-fluid/CDM coefficient g=-43/64, closure=0 and legacy identity=0.
  const auto legacy = prepare_thermal_background(
      {70, 3. / 32, 1. / 8, 0, 17. / 32, {}});
  const auto legacy_epoch = detail::thermal_conformal_epoch(
      legacy, .5L, .001L, 7.L / 32, policy);
  need(legacy_epoch.status == S::ok && legacy_epoch.fb == 0 &&
           legacy_epoch.fc == 17.L / 32 && legacy_epoch.g == -43.L / 64 &&
           legacy_epoch.closure_defect == 0 &&
           legacy_epoch.background_identity_defect == 0 &&
           legacy_epoch.hcal == epoch.hcal && legacy_epoch.x2 == epoch.x2,
       "original zero-baryon fluid arithmetic and retained H are preserved");
  need(prepare_perfect_fluid_transfer(legacy, 1e-9).status() == S::ok,
       "original zero-baryon public domain remains available");
}
void actual_failed_query() {
  // Positive-mass thermal source is used ONLY to obtain a real charged
  // momentum failure. This does not admit its perturbations in any consumer.
  const auto background = prepare_thermal_background(
      {70, .0001, 0, .04, .2, {{.06, .00017, 2}}});
  need(background.status() == S::ok, "failed-query thermal owner");
  ThermalPolicy capped;
  capped.momentum_method = background.momentum_method();
  capped.maximum_callbacks_per_evaluation = 5;
  capped.maximum_total_callbacks = 5;
  const auto failed = background.scaled_expansion(.5L, capped);
  need(failed.status != S::ok && failed.callbacks > 0 &&
           failed.callbacks <= capped.maximum_total_callbacks,
       "actual bounded momentum failure retains attempted work");
  const auto imported = detail::thermal_conformal_epoch_from_scaled(
      background, failed, .5L, .001L, .0001L);
  need(imported.status == failed.status, "actual failed import status retained");
  raw_matches(imported, failed);
  const auto wrapped = detail::thermal_conformal_epoch(
      background, .5L, .001L, .0001L, capped);
  need(wrapped.status == failed.status, "actual failed wrapper status retained");
  raw_matches(wrapped, failed);
  const auto preflight = detail::thermal_conformal_epoch(
      background, 0, .001L, .0001L, capped);
  need(preflight.status == S::invalid_input && !preflight.raw_scaled_query &&
           preflight.momentum_callbacks == 0 && !preflight.forward &&
           !preflight.shadow,
       "wrapper preflight refusal records no attempted query");
  std::cout << std::setprecision(std::numeric_limits<W>::max_digits10)
            << "actual_failed_query status=" << static_cast<unsigned>(failed.status)
            << " callbacks=" << failed.callbacks << " p=" << failed.a4_e2
            << " p_error=" << failed.error_estimate << '\n';
}
void retained_coefficient_capture() {
  need(!detail::ThermalRetainedCoefficientAccess::capture(ThermalBackground{}),
       "default background has no retained coefficient witness");
  const auto failed = prepare_thermal_background({70, 0, 0, 0, 2, {}});
  need(failed.status() == S::outside_domain &&
           !detail::ThermalRetainedCoefficientAccess::capture(failed),
       "failed physical closure has no retained coefficient witness");
  // Binary64 .3 is below 3/10. The exact retained 1-Omega_c rounds down to
  // emitted binary64 .7, with the exact positive private-minus-getter loss2^-54.
  const auto background = prepare_thermal_background({70, 0, 0, 0, .3, {}});
  const auto witness = detail::ThermalRetainedCoefficientAccess::capture(background);
  need(witness && witness->status == S::ok &&
           witness->lambda_retained == 1 - W(.3) &&
           witness->lambda_emitted == .7 &&
           witness->lambda_getter_signed_loss == 0x1p-54L &&
           witness->omega_species_today == 0 &&
           witness->species_normalization_error == 0 &&
           witness->critical_density_ev4 > 0 &&
           std::isfinite(witness->critical_density_ev4) &&
           witness->momentum_method == background.momentum_method(),
       "actual retained coefficient and exact getter cast-loss witness");
  const auto endpoint = detail::thermal_conformal_epoch(
      background, 1, .001L, 0, ThermalPolicy{});
  need(endpoint.status == S::ok && endpoint.p == 1 && endpoint.p_error == 0 &&
           endpoint.closure_defect == -0x1p-54L && !endpoint.forward &&
           !endpoint.shadow,
       "normalized P endpoint does not reclose public getter fractions");
  auto copied = background;
  auto moved = std::move(copied);
  need(!detail::ThermalRetainedCoefficientAccess::capture(copied),
       "moved background has no coefficient witness");
  const auto moved_witness = detail::ThermalRetainedCoefficientAccess::capture(moved);
  need(moved_witness && moved_witness->lambda_retained == witness->lambda_retained &&
           moved_witness->lambda_getter_signed_loss ==
               witness->lambda_getter_signed_loss,
       "copied then moved background retains matching scalar state");
  const auto rounding = std::fegetround();
  need(rounding == FE_TONEAREST && std::fesetround(FE_UPWARD) == 0,
       "capture rounding adversary available");
  const bool refused = !detail::ThermalRetainedCoefficientAccess::capture(background);
  need(std::fesetround(rounding) == 0 && refused,
       "capture checks current arithmetic before getter conversion");
}
void no_species_shadow_diagnostic() {
  const auto mapping = map_thermal_physical_model({70, .02, .10, 2.7, 0, {}});
  need(mapping.status == S::ok && mapping.model && mapping.scalar_witnesses,
       "actual one physical map retains all scalar provenance");
  const auto background = prepare_thermal_background(*mapping.model);
  const auto captured = detail::ThermalRetainedCoefficientAccess::capture(background);
  need(captured.has_value(), "actual same retained coefficient capture");
  for (W a : {1e-10L, .01L, 1.L}) {
    auto epoch = detail::thermal_conformal_epoch(background, a, .01L,
        mapping.model->omega_gamma, ThermalPolicy{});
    const auto p = epoch.p, h = epoch.hcal, g = epoch.g;
    const auto raw = epoch.raw_scaled_query;
    need(detail::thermal_conformal_no_species_diagnostics(background, *captured,
             *mapping.scalar_witnesses, a, epoch) == S::ok && epoch.forward && epoch.shadow,
         "conditional no-species forward/shadow provenance available");
    const auto &s = *epoch.shadow; const auto &f = *epoch.forward;
    need(epoch.p == p && epoch.hcal == h && epoch.g == g && raw &&
         epoch.raw_scaled_query->a4_e2 == raw->a4_e2 &&
         epoch.raw_scaled_query->callbacks == raw->callbacks,
         "diagnostics reuse actual P/H and raw query without substitution");
    need(s.p_shadow > s.p_shadow_radius && s.actual_p_source_radius > 0 &&
         s.q == s.p_shadow_n / s.p_shadow && s.ell == s.q / 2 - 1 &&
         s.dp == s.q - 2 * (g + 1) && s.q_radius > 0 && s.dp_radius > 0 &&
         std::abs(g - s.ell) <= f.g_absolute_estimate &&
         f.hcal_absolute_estimate > 0 && f.hcal_absolute_estimate < h &&
         f.fr_absolute_estimate == f.fg_absolute_estimate,
         "literal shadow derivative and actual coefficient discrepancy retained");
    if (a == 1) need(epoch.p == 1 && epoch.p_error == 0,
                     "original P1 normalized shortcut remains literal");
    auto mismatch = *mapping.scalar_witnesses;
    mismatch[1].emitted_value = .5;
    need(detail::thermal_conformal_no_species_diagnostics(background, *captured,
             mismatch, a, epoch) == S::conditioning_budget_exceeded &&
         !epoch.forward && !epoch.shadow && epoch.p == p && epoch.hcal == h,
         "unmatched map clears dependencies without erasing actual state");
  }
  auto absent = detail::ThermalConformalEpoch{};
  need(detail::thermal_conformal_no_species_diagnostics(background, *captured,
           *mapping.scalar_witnesses, .5L, absent) == S::conditioning_budget_exceeded &&
       !absent.forward && !absent.shadow, "missing charged query cannot become a zero envelope");
  auto zero_k = detail::thermal_conformal_epoch(background, .5L, 0,
      mapping.model->omega_gamma, ThermalPolicy{});
  need(zero_k.status == S::ok && zero_k.x2 == 0 &&
       detail::thermal_conformal_no_species_diagnostics(background, *captured,
           *mapping.scalar_witnesses, .5L, zero_k) == S::outside_domain &&
       !zero_k.forward && !zero_k.shadow,
       "positive-k diagnostic slice refuses zero x2 while old conversion remains available");
  const auto rounding = std::fegetround();
  auto epoch = detail::thermal_conformal_epoch(background, .5L, .01L,
      mapping.model->omega_gamma, ThermalPolicy{});
  need(std::fesetround(FE_UPWARD) == 0, "diagnostic rounding adversary available");
  const auto refused = detail::thermal_conformal_no_species_diagnostics(
      background, *captured, *mapping.scalar_witnesses, .5L, epoch);
  need(std::fesetround(rounding) == 0 && refused == S::conditioning_budget_exceeded &&
       !epoch.forward && !epoch.shadow, "current nonnearest diagnostic refuses");
}
} // namespace
int main() {
  dyadic_source();
  actual_failed_query();
  retained_coefficient_capture();
  no_species_shadow_diagnostic();
  std::cout << "thermal_conformal_epoch_contract passed "
               "dyadic/totalmatter/legacy/query-receipts/absence controls "
            << "epoch_bytes=" << sizeof(irred::cosmology::detail::ThermalConformalEpoch)
            << '\n';
}
