#pragma once
#include "irred/thermal_neutrino.hpp"
#include "flat_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace irred::cosmology::detail {
using ThermalWide = long double;
inline constexpr ThermalWide thermal_loading_floor = 64 * std::numeric_limits<ThermalWide>::epsilon();
struct ThermalBaryonLoading {
  numerics::Status status = numerics::Status::invalid_input;
  ThermalWide ratio_today = 0, arithmetic_estimate = 0;
  ThermalWide numerator = 0, denominator = 0;
};
inline ThermalBaryonLoading thermal_baryon_loading(const ThermalBackground &background) noexcept {
  using S = numerics::Status;
  ThermalBaryonLoading out;
  if (background.status() != S::ok) { out.status = background.status(); return out; }
  const auto &m = background.source();
  if (!(m.omega_gamma > 0) || !(m.omega_b >= 0)) { out.status = S::outside_domain; return out; }
  out.numerator = 3 * ThermalWide(m.omega_b);
  out.denominator = 4 * ThermalWide(m.omega_gamma);
  out.ratio_today = out.numerator / out.denominator;
  // Frozen ratio law: 32 wide eps plus the measured binary64 reporting cast.
  // Additional helper operations have their separate 64-eps diagnostics below.
  out.arithmetic_estimate = 32 * std::numeric_limits<ThermalWide>::epsilon() * std::abs(out.ratio_today) +
      std::abs(out.ratio_today - ThermalWide(static_cast<double>(out.ratio_today)));
  if (!std::isfinite(out.ratio_today) || !std::isfinite(out.arithmetic_estimate) ||
      (out.ratio_today != 0 && (!std::isnormal(out.ratio_today) || !std::isnormal(out.arithmetic_estimate)))) {
    out.status = S::outside_domain; return out;
  }
  out.status = S::ok; return out;
}
struct ThermalLoadingValue {
  numerics::Status status = numerics::Status::invalid_input;
  ThermalWide value = 0, arithmetic_estimate = 0;
};
inline ThermalLoadingValue thermal_baryon_loading_at_inverse_scale(
    const ThermalBaryonLoading &loading, ThermalWide u) noexcept {
  ThermalLoadingValue out;
  if (loading.status != numerics::Status::ok) { out.status = loading.status; return out; }
  if (!(u > 0) || !std::isfinite(u)) return out;
  if (loading.numerator == 0) { out.status = numerics::Status::ok; return out; }
  // Preserve the original pure-H operation order exactly. A represented zero
  // from positive numerator is a refusal, even if the finite u is enormous.
  out.value = loading.numerator / (loading.denominator * u);
  out.arithmetic_estimate = thermal_loading_floor * out.value;
  out.status = (out.value > 0 && std::isnormal(out.value) && std::isnormal(out.arithmetic_estimate))
      ? numerics::Status::ok : numerics::Status::outside_domain;
  return out;
}
struct ThermalBaryonSound {
  numerics::Status status = numerics::Status::invalid_input;
  ThermalWide loading = 0, loading_estimate = 0, denominator = 0;
  ThermalWide sound_speed_over_c = 0, sound_speed_estimate_over_c = 0;
  ThermalWide sound_speed_squared_over_c_squared = 0;
  ThermalWide sound_speed_squared_estimate_over_c_squared = 0;
};
inline ThermalBaryonSound thermal_baryon_sound(const ThermalBaryonLoading &r, ThermalWide a) noexcept {
  using S = numerics::Status;
  ThermalBaryonSound out;
  if (r.status != S::ok) { out.status = r.status; return out; }
  // Shared early/late consumers own only 0<=a<=1. Future scales need a
  // separately reviewed domain, rather than overflow-shaped zero sound speed.
  if (!std::isfinite(a)) { out.status = S::nonfinite_input; return out; }
  if (!(a >= 0) || a > 1) { out.status = S::outside_domain; return out; }
  out.loading = r.ratio_today * a;
  out.loading_estimate = r.arithmetic_estimate * a + thermal_loading_floor * out.loading;
  out.denominator = 1 + out.loading;
  const ThermalWide e = out.loading_estimate + thermal_loading_floor * out.denominator;
  if (!(out.denominator > e) || !std::isfinite(out.denominator)) { out.status = S::conditioning_budget_exceeded; return out; }
  out.sound_speed_squared_over_c_squared = 1 / (3 * out.denominator);
  out.sound_speed_over_c = 1 / std::sqrt(3 * out.denominator);
  const ThermalWide lo = out.denominator - e, hi = out.denominator + e;
  out.sound_speed_squared_estimate_over_c_squared = std::max(
      out.sound_speed_squared_over_c_squared - 1 / (3 * hi),
      1 / (3 * lo) - out.sound_speed_squared_over_c_squared) + thermal_loading_floor * out.sound_speed_squared_over_c_squared;
  // Stable positive inverse-square-root difference; no cancellation at tiny e.
  out.sound_speed_estimate_over_c = out.sound_speed_over_c *
      e / (std::sqrt(lo) * (std::sqrt(out.denominator) + std::sqrt(lo))) +
      thermal_loading_floor * out.sound_speed_over_c;
  if (!std::isnormal(out.sound_speed_over_c) || !(out.sound_speed_over_c > 0) ||
      !std::isnormal(out.sound_speed_squared_over_c_squared) || !(out.sound_speed_squared_over_c_squared > 0) ||
      !std::isnormal(out.sound_speed_estimate_over_c) || !(out.sound_speed_estimate_over_c > 0) ||
      !std::isnormal(out.sound_speed_squared_estimate_over_c_squared) || !(out.sound_speed_squared_estimate_over_c_squared > 0)) {
    out.status = S::conditioning_budget_exceeded; return out;
  }
  out.status = S::ok; return out;
}

struct ThermalRulerAllowance {
  ThermalPolicy thermal;
  double absolute_tolerance_mpc = 0, relative_tolerance = 0;
  std::size_t maximum_outer_callbacks = 0;
  unsigned maximum_depth = 30;
};
// Consumer-owned spent work and category fields are checked before every add.
// A new scoped integral never resets aggregate/history work. The constructor
// admits only bounded counters; public output records cannot supply freshness.
class ThermalRulerWorkBudget {
  std::size_t &spent_, &outer_, &momentum_;
  std::size_t cap_, start_, per_call_, momentum_cap_;
  bool valid_;
public:
  ThermalRulerWorkBudget(std::size_t &spent, std::size_t cap, std::size_t &outer,
                        std::size_t &momentum, std::size_t per_call,
                        std::size_t momentum_cap) noexcept
      : spent_(spent), outer_(outer), momentum_(momentum), cap_(cap), start_(spent),
        per_call_(per_call), momentum_cap_(momentum_cap),
        valid_(spent <= cap && momentum <= momentum_cap && outer <= spent && momentum <= spent) {}
  std::size_t remaining() const noexcept {
    if (!valid_ || spent_ > cap_ || spent_ < start_ || spent_ - start_ > per_call_) return 0;
    return std::min(cap_ - spent_, per_call_ - (spent_ - start_));
  }
  std::size_t momentum_remaining() const noexcept {
    return valid_ && momentum_ <= momentum_cap_ ? std::min(remaining(), momentum_cap_ - momentum_) : 0;
  }
  bool charge_outer() noexcept {
    if (!remaining() || outer_ == SIZE_MAX) return false;
    ++spent_; ++outer_; return true;
  }
  bool charge_momentum(std::size_t n) noexcept {
    if (n > momentum_remaining()) { valid_ = false; return false; }
    spent_ += n; momentum_ += n; return true;
  }
};
struct ThermalRulerIntegral {
  numerics::Status status = numerics::Status::invalid_input;
  ThermalWide value_mpc = 0, error_estimate_mpc = 0;
  ThermalWide outer_quadrature_estimate_mpc = 0, background_estimate_mpc = 0;
  ThermalWide loading_estimate_mpc = 0, arithmetic_estimate_mpc = 0;
};
namespace thermal_ruler_internal {
struct Context {
  const ThermalBackground &background;
  const ThermalBaryonLoading &loading;
  const ThermalRulerAllowance &allowance;
  ThermalRulerWorkBudget &work;
  ThermalWide endpoint, background_relative = 0, kernel_cast_relative = 0;
  std::size_t outer_attempts = 0;
  numerics::Status status = numerics::Status::ok;
};
inline double integrand(double t, const void *ptr) {
  auto &c = *const_cast<Context *>(static_cast<const Context *>(ptr));
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  if (c.status != numerics::Status::ok) return nan;
  if (c.outer_attempts >= c.allowance.maximum_outer_callbacks || !c.work.charge_outer()) {
    c.status = numerics::Status::work_limit; return nan;
  }
  ++c.outer_attempts;
  auto p = c.allowance.thermal;
  p.maximum_total_callbacks = std::min(p.maximum_total_callbacks, c.work.momentum_remaining());
  const ThermalWide a = c.endpoint * t;
  const auto q = c.background.scaled_expansion(a, p);
  if (!c.work.charge_momentum(q.callbacks)) { c.status = numerics::Status::work_limit; return nan; }
  c.status = q.status;
  if (c.status != numerics::Status::ok) return nan;
  if (!(q.a4_e2 > q.error_estimate) || !(q.error_estimate >= 0)) {
    c.status = numerics::Status::conditioning_budget_exceeded; return nan;
  }
  const ThermalWide low = q.a4_e2 - q.error_estimate;
  const ThermalWide relative = q.error_estimate /
      (std::sqrt(low) * (std::sqrt(q.a4_e2) + std::sqrt(low)));
  c.background_relative = std::max(c.background_relative, relative + thermal_loading_floor * (1 + relative));
  const auto sound = thermal_baryon_sound(c.loading, a);
  if (sound.status != numerics::Status::ok) { c.status = sound.status; return nan; }
  // Original ruler center and operation order; cs itself is not substituted.
  const ThermalWide value = 1 / std::sqrt(q.a4_e2 * sound.denominator);
  const double rounded = static_cast<double>(value);
  if (!(rounded > 0) || !physical_representable(value)) { c.status = numerics::Status::outside_domain; return nan; }
  c.kernel_cast_relative = std::max(c.kernel_cast_relative, std::abs(value - rounded) / value);
  return rounded;
}
// Root-pinned stable loading expression. Arithmetic estimates are diagnostics,
// not a rigorous libm or continuous physical enclosure.
inline std::optional<ThermalWide> loading_xi(const ThermalBaryonLoading &r, ThermalWide a) noexcept {
  if (r.arithmetic_estimate == 0 || a == 0) return 0;
  const ThermalWide A = 1 + r.ratio_today * a,
      B = 1 + (r.ratio_today - r.arithmetic_estimate) * a,
      product = r.arithmetic_estimate * a;
  if (!(B > 0) || !std::isnormal(product) || !std::isfinite(A)) return {};
  const ThermalWide xi = product / (std::sqrt(B) * (std::sqrt(A) + std::sqrt(B)));
  if (!std::isnormal(xi)) return {};
  return xi + thermal_loading_floor * xi;
}
} // namespace thermal_ruler_internal
inline ThermalRulerIntegral integrate_thermal_ruler(
    const ThermalBackground &background, const ThermalBaryonLoading &loading,
    ThermalWide endpoint, const ThermalRulerAllowance &p, ThermalRulerWorkBudget &work) {
  using S = numerics::Status;
  ThermalRulerIntegral out;
  if (loading.status != S::ok) { out.status = loading.status; return out; }
  if (!(endpoint > 0) || endpoint > 1 || !std::isfinite(endpoint) ||
      p.maximum_depth > 60 || p.absolute_tolerance_mpc < 0 || p.relative_tolerance < 0 ||
      !std::isfinite(p.absolute_tolerance_mpc) || !std::isfinite(p.relative_tolerance)) return out;
  const auto xi = thermal_ruler_internal::loading_xi(loading, endpoint);
  if (!xi) { out.status = S::conditioning_budget_exceeded; return out; }
  const ThermalWide scale = prepare_flat_scale(background.source().h0_km_s_mpc).distance_mpc * endpoint / std::sqrt(3.L);
  const double absolute = static_cast<double>(std::min(ThermalWide(p.absolute_tolerance_mpc) / (8 * scale), ThermalWide(std::numeric_limits<double>::max())));
  if (!(absolute > 0 || p.relative_tolerance > 0)) { out.status = S::conditioning_budget_exceeded; return out; }
  const auto available = std::min(work.remaining(), p.maximum_outer_callbacks);
  if (available < 3) { out.status = S::work_limit; return out; }
  thermal_ruler_internal::Context c{background, loading, p, work, endpoint};
  const auto q = numerics::integrate(thermal_ruler_internal::integrand, &c, 0, 1,
      {absolute, p.relative_tolerance / 8, available, p.maximum_depth});
  out.status = c.status == S::ok ? q.status : c.status;
  if (out.status != S::ok) return out;
  out.value_mpc = scale * q.value;
  out.outer_quadrature_estimate_mpc = scale * q.error_estimate;
  if (!(c.kernel_cast_relative < 1)) { out.status = S::conditioning_budget_exceeded; return out; }
  const ThermalWide response_scale = (std::abs(out.value_mpc) + out.outer_quadrature_estimate_mpc) /
      (1 - c.kernel_cast_relative);
  out.background_estimate_mpc = response_scale * c.background_relative;
  const ThermalWide cast_component = response_scale * c.kernel_cast_relative,
      operation_component = out.value_mpc * (64 * std::numeric_limits<double>::epsilon());
  out.arithmetic_estimate_mpc = cast_component + operation_component +
      scale * std::abs(ThermalWide(std::nextafter(q.value, INFINITY)) - q.value) / 2;
  const ThermalWide E0 = out.outer_quadrature_estimate_mpc + out.background_estimate_mpc + out.arithmetic_estimate_mpc;
  out.loading_estimate_mpc = (std::abs(out.value_mpc) + E0) * *xi;
  out.error_estimate_mpc = E0 + out.loading_estimate_mpc;
  const bool positive_components_resolved =
      (q.error_estimate == 0 || std::isnormal(out.outer_quadrature_estimate_mpc)) &&
      (c.background_relative == 0 || std::isnormal(out.background_estimate_mpc)) &&
      (c.kernel_cast_relative == 0 || std::isnormal(cast_component)) &&
      std::isnormal(operation_component) &&
      (*xi == 0 || std::isnormal(out.loading_estimate_mpc));
  if (!positive_components_resolved || !std::isfinite(out.error_estimate_mpc) || !std::isnormal(out.value_mpc) ||
      !(out.error_estimate_mpc > 0) || !std::isnormal(out.error_estimate_mpc)) out.status = S::conditioning_budget_exceeded;
  return out;
}
struct ThermalSoundDistanceBound {
  numerics::Status status = numerics::Status::invalid_input;
  ThermalWide upper_mpc_per_redshift = 0;
};
inline ThermalSoundDistanceBound thermal_sound_distance_bound(
    const ThermalBackground &background, const ThermalBaryonLoading &r,
    ThermalWide amin, ThermalWide amax, const ThermalPolicy &policy,
    ThermalRulerWorkBudget &work) {
  using S = numerics::Status;
  ThermalSoundDistanceBound out;
  if (r.status != S::ok) { out.status = r.status; return out; }
  if (!(amin > 0) || !(amax >= amin) || amax > 1 || !background.source().species.empty()) return out;
  const auto lambda = background.omega_lambda();
  if (!lambda || *lambda < 0) { out.status = S::outside_domain; return out; }
  if (!work.charge_outer()) { out.status = S::work_limit; return out; }
  auto p = policy; p.maximum_total_callbacks = std::min(p.maximum_total_callbacks, work.momentum_remaining());
  const auto q = background.scaled_expansion(amin, p);
  if (!work.charge_momentum(q.callbacks)) { out.status = S::work_limit; return out; }
  out.status = q.status;
  if (out.status != S::ok) return out;
  const ThermalWide lowerP = q.a4_e2 - q.error_estimate,
      lowerR = 1 + (r.ratio_today - r.arithmetic_estimate) * amin;
  if (!(lowerP > 0) || !(lowerR > 0)) { out.status = S::conditioning_budget_exceeded; return out; }
  out.upper_mpc_per_redshift = prepare_flat_scale(background.source().h0_km_s_mpc).distance_mpc *
      amax * amax / std::sqrt(3 * lowerP * lowerR);
  out.upper_mpc_per_redshift *= 1 + thermal_loading_floor;
  if (!std::isnormal(out.upper_mpc_per_redshift)) out.status = S::conditioning_budget_exceeded;
  return out;
}
} // namespace irred::cosmology::detail
