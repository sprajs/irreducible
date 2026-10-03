#include "irred/effective_fluid.hpp"
#include "payload_accounting.hpp"
#include "thermal_conformal_epoch.hpp"
#include "thermal_ruler.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <utility>
namespace irred::cosmology {
namespace {
using S = numerics::Status;
using W = long double;
constexpr W floor = 128 * std::numeric_limits<W>::epsilon();
bool profile() noexcept {
  return std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384 && std::fegetround() == FE_TONEAREST;
}
bool normal(W x) noexcept { return std::isfinite(x) && (x == 0 || std::isnormal(x)); }
bool tolerance(double a, double r) noexcept {
  return std::isfinite(a) && std::isfinite(r) && a >= 0 && r >= 0 && (a > 0 || r > 0);
}
bool valid(const ThermalObservablePolicy &p) noexcept {
  return profile() && tolerance(p.absolute_tolerance_mpc, p.relative_tolerance) &&
         tolerance(p.thermal.absolute_tolerance, p.thermal.relative_tolerance) &&
         p.maximum_depth <= 60 && p.thermal.maximum_depth <= 60 &&
         !thermal_momentum_method_id(p.thermal.momentum_method).empty();
}
std::optional<std::size_t> origins(const EffectiveFluidRequest &r) noexcept {
  std::size_t n = 0;
  if (!irred::detail::checked_payload_add(n, r.source_origin.size(), 1) ||
      !irred::detail::checked_payload_add(n, r.endpoint_origin.size(), 1) ||
      !irred::detail::checked_payload_add(n, 2, 1)) return {};
  return n;
}
W power(W x) noexcept { const W x2 = x * x; return x2 * x2 * std::sqrt(x); }
EffectiveFluidState shape(W xc, W ac, W a) noexcept {
  EffectiveFluidState out;
  if (!profile()) return out;
  if (!std::isfinite(a)) { out.status = S::nonfinite_input; return out; }
  if (a < 0 || a > 1) { out.status = S::outside_domain; return out; }
  W fraction = 0;
  if (a != 0) {
    const W ratio = a <= ac ? a / ac : ac / a;
    const W p = power(ratio);
    if (!(p > 0) || !normal(p)) { out.status = S::outside_domain; return out; }
    fraction = a <= ac ? p / (1 + p) : 1 / (1 + p);
    out.density = xc == 0 ? 0 : (a <= ac ? 2 * xc / (1 + p) : 2 * xc * p / (1 + p));
  } else out.density = 2 * xc;
  out.equation_of_state = -1 + 1.5L * fraction;
  out.log_scale_derivative = -4.5L * out.density * fraction;
  out.density_estimate = floor * out.density;
  out.equation_of_state_estimate = floor * (1 + std::abs(out.equation_of_state));
  out.derivative_estimate = floor * std::abs(out.log_scale_derivative);
  if (!normal(out.density) || !normal(out.density_estimate) ||
      !normal(out.log_scale_derivative) || !normal(out.derivative_estimate) ||
      (xc > 0 && (!(out.density > 0) || !(out.density_estimate > 0))) ||
      (xc > 0 && a > 0 && (!(out.log_scale_derivative < 0) || !(out.derivative_estimate > 0)))) {
    out.status = S::outside_domain; return out;
  }
  out.status = S::ok; return out;
}
void project(ThermalBackgroundValue &out, W value, W error, double relative,
             bool exact) {
  if (!(value > 0) || !detail::physical_representable(value)) { out.status = S::outside_domain; return; }
  const double rounded = static_cast<double>(value);
  const W diagnostic = exact ? 0 : error + 64 * std::numeric_limits<double>::epsilon() * value + std::abs(value - rounded);
  if (!detail::physical_representable(diagnostic) || diagnostic > relative * value) {
    out.status = S::conditioning_budget_exceeded; return;
  }
  out = {S::ok, rounded, static_cast<double>(diagnostic)};
}
} // namespace
std::optional<std::size_t> effective_fluid_payload_bound(
    std::size_t points, std::size_t species, std::size_t origin_bytes) noexcept {
  if (species > 16 || points > 65536) return {};
  const auto thermal = thermal_background_payload_bound(points, species);
  if (!thermal) return {};
  irred::detail::PayloadAccounting b(*thermal);
  b.add(1, 2 * sizeof(EffectiveFluidBackground) + sizeof(EffectiveFluidBatch) +
      sizeof(detail::ThermalRetainedCoefficientWitness) + 4096);
  b.add(species, 4 * sizeof(ThermalSpecies));
  b.add(points, 2 * sizeof(ThermalBackgroundRow));
  b.add(origin_bytes, 4);
  return b.result();
}
EffectiveFluidBackground prepare_effective_fluid(const EffectiveFluidRequest &r,
                                                ThermalObservablePolicy p) {
  EffectiveFluidBackground out;
  if (!valid(p)) return out;
  const auto n = origins(r);
  const auto bytes = n ? effective_fluid_payload_bound(0, r.reference.species.size(), *n) : std::nullopt;
  if (!bytes) return out;
  if (*bytes > p.maximum_native_bytes || r.reference.species.size() > p.thermal.maximum_species) {
    out.status_ = S::work_limit; return out;
  }
  if (!std::isfinite(r.density_at_transition) || !std::isfinite(r.transition_scale) ||
      !std::isfinite(r.ruler_endpoint)) { out.status_ = S::nonfinite_input; return out; }
  if (r.density_at_transition < 0 || !(r.transition_scale > 0 && r.transition_scale < 1) ||
      !(r.ruler_endpoint > 0 && r.ruler_endpoint <= 1) ||
      r.source_origin.empty() || r.endpoint_origin.empty()) { out.status_ = S::outside_domain; return out; }
  out.request_ = r; // Allocate all consumer source storage before FD work.
  auto nested = p.thermal;
  nested.maximum_total_callbacks = std::min(nested.maximum_total_callbacks, p.maximum_total_callbacks);
  out.reference_ = prepare_thermal_background(r.reference, nested);
  out.status_ = out.reference_.status();
  if (out.status_ != S::ok) return out;
  const auto retained = detail::ThermalRetainedCoefficientAccess::capture(out.reference_);
  if (!retained) { out.status_ = S::conditioning_budget_exceeded; return out; }
  if (!(r.reference.omega_gamma > 0)) { out.status_ = S::outside_domain; return out; }
  out.lambda_ = retained->lambda_retained;
  if (r.density_at_transition == 0) return out;
  const auto today = shape(r.density_at_transition, r.transition_scale, 1);
  out.status_ = today.status;
  if (out.status_ != S::ok) return out;
  out.today_density_ = today.density; out.today_estimate_ = today.density_estimate;
  out.transition_power_ = power(r.transition_scale);
  if (!(out.transition_power_ > 0) || !normal(out.transition_power_) || out.today_density_ > out.lambda_) {
    out.status_ = S::outside_domain; return out;
  }
  out.lambda_ -= out.today_density_;
  const W closure_estimate = retained->species_normalization_error + out.today_estimate_ +
      floor * (retained->lambda_retained + out.today_density_);
  if (!normal(out.lambda_) || !normal(closure_estimate) || out.lambda_ < closure_estimate)
    out.status_ = S::conditioning_budget_exceeded;
  return out;
}
EffectiveFluidBackground::EffectiveFluidBackground(EffectiveFluidBackground &&other) noexcept { *this = std::move(other); }
EffectiveFluidBackground &EffectiveFluidBackground::operator=(const EffectiveFluidBackground &other) {
  if (this != &other) {
    // Replace old capacities together with the complete new physical owner.
    auto replacement = other;
    *this = std::move(replacement);
  }
  return *this;
}
EffectiveFluidBackground &EffectiveFluidBackground::operator=(EffectiveFluidBackground &&other) noexcept {
  if (this == &other) return *this;
  status_ = other.status_; request_ = std::move(other.request_); reference_ = std::move(other.reference_);
  today_density_ = other.today_density_; today_estimate_ = other.today_estimate_;
  transition_power_ = other.transition_power_; lambda_ = other.lambda_;
  other.status_ = S::invalid_input; other.request_ = {};
  other.today_density_ = other.today_estimate_ = other.transition_power_ = other.lambda_ = 0;
  return *this;
}
std::optional<long double> EffectiveFluidBackground::omega_lambda() const noexcept {
  if (status_ != S::ok || !profile()) return {};
  return lambda_;
}
EffectiveFluidState EffectiveFluidBackground::fluid_state(W a) const noexcept {
  if (status_ != S::ok) { EffectiveFluidState out; out.status = status_; return out; }
  return shape(request_.density_at_transition, request_.transition_scale, a);
}
ThermalScaledExpansion EffectiveFluidBackground::scaled_expansion(W a, ThermalPolicy p) const {
  ThermalScaledExpansion out;
  if (status_ != S::ok) { out.status = status_; return out; }
  out = reference_.scaled_expansion(a, p);
  if (out.status != S::ok || request_.density_at_transition == 0 || a == 0 || a == 1) return out;
  const auto fluid = fluid_state(a);
  if (fluid.status != S::ok) { out.status = fluid.status; return out; }
  const W factor = -std::expm1(4.5L * std::log(a)),
      correction = fluid.density * factor / (1 + transition_power_),
      a2 = a * a, a4 = a2 * a2, addition = a4 * correction;
  const W addition_error = a4 * (fluid.density_estimate * factor / (1 + transition_power_) + floor * correction);
  if (!(factor > 0) || !(correction > 0) || !(addition > 0) || !(addition_error > 0) ||
      !normal(factor) || !normal(correction) || !normal(addition) || !normal(addition_error)) {
    out.status = S::outside_domain; return out;
  }
  out.a4_e2 += addition;
  out.error_estimate += addition_error + floor * out.a4_e2;
  if (!normal(out.a4_e2) || !normal(out.error_estimate) || !(out.a4_e2 > out.error_estimate))
    out.status = S::conditioning_budget_exceeded;
  return out;
}
EffectiveFluidBatch EffectiveFluidBackground::evaluate(std::span<const double> factors,
    unsigned mask, ThermalObservablePolicy p) const {
  EffectiveFluidBatch out;
  if (!valid(p) || !mask || (mask & ~(thermal_e | thermal_h | effective_fluid_ruler))) return out;
  if (status_ != S::ok) { out.status = status_; return out; }
  if (p.thermal.momentum_method != reference_.momentum_method()) return out;
  const auto n = origins(request_);
  const auto bytes = n ? effective_fluid_payload_bound(factors.size(), source().species.size(), *n) : std::nullopt;
  if (!bytes) return out;
  if (*bytes > p.maximum_native_bytes || factors.size() > p.maximum_points ||
      source().species.size() > p.thermal.maximum_species) { out.status = S::work_limit; return out; }
  out.rows.reserve(factors.size()); out.requested_outputs = mask; out.status = S::ok;
  if (mask & effective_fluid_ruler) {
    out.ruler = EarlyLateValue{};
    const auto loading = detail::thermal_baryon_loading(reference_);
    detail::ThermalRulerWorkBudget work(out.callbacks, p.maximum_total_callbacks,
        out.outer_callbacks, out.momentum_callbacks, p.maximum_callbacks_per_point,
        p.thermal.maximum_total_callbacks);
    const detail::ThermalRulerAllowance allowance{p.thermal, p.absolute_tolerance_mpc,
        p.relative_tolerance, p.maximum_callbacks_per_point, p.maximum_depth};
    const auto q = detail::integrate_thermal_ruler(*this, loading, W(request_.ruler_endpoint), allowance, work);
    out.ruler->status = q.status;
    if (q.status == S::ok) {
      const double rounded = static_cast<double>(q.value_mpc);
      const W error = q.error_estimate_mpc + 64 * std::numeric_limits<double>::epsilon() * q.value_mpc + std::abs(q.value_mpc - rounded);
      if (!detail::physical_representable(q.value_mpc) || !detail::physical_representable(error)) out.ruler->status = S::outside_domain;
      else if (error > p.absolute_tolerance_mpc + p.relative_tolerance * q.value_mpc) out.ruler->status = S::conditioning_budget_exceeded;
      else *out.ruler = {S::ok, rounded, static_cast<double>(error)};
    }
  }
  for (double a : factors) {
    out.rows.emplace_back(); auto &row = out.rows.back(); row.scale_factor = a;
    if (!(mask & (thermal_e | thermal_h))) continue;
    S status = S::ok;
    if (!std::isfinite(a)) status = S::nonfinite_input;
    else if (!(a > 0 && a <= 1)) status = S::outside_domain;
    W e = 0, error = 0;
    if (status == S::ok) {
      auto nested = p.thermal;
      nested.maximum_total_callbacks = std::min({p.maximum_callbacks_per_point,
          p.maximum_total_callbacks - out.callbacks, p.thermal.maximum_total_callbacks - out.momentum_callbacks});
      const auto q = scaled_expansion(a, nested);
      row.callbacks = q.callbacks; out.momentum_callbacks += q.callbacks; out.callbacks += q.callbacks;
      status = q.status;
      if (status == S::ok) {
        const W a2 = W(a) * a;
        e = std::sqrt(q.a4_e2) / a2;
        error = q.error_estimate / (2 * std::sqrt(q.a4_e2 - q.error_estimate) * a2);
      }
    }
    if (mask & thermal_e) {
      row.e.status = status;
      if (status == S::ok) project(row.e, e, error, p.relative_tolerance, a == 1);
    }
    if (mask & thermal_h) {
      row.h_km_s_mpc.status = status;
      if (status == S::ok) project(row.h_km_s_mpc, e * source().h0_km_s_mpc,
          error * source().h0_km_s_mpc, p.relative_tolerance, a == 1);
    }
  }
  return out;
}
} // namespace irred::cosmology
