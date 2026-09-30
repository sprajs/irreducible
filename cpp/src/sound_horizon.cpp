#include "irred/sound_horizon.hpp"
#include "flat_geometry.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <limits>

namespace irred::cosmology {
namespace {
struct State {
  long double radiation, matter, lambda, baryon_loading, a_drag;
};
numerics::Status state(const SoundHorizonRequest &r, State &s) {
  const auto &m = r.model;
  const double values[]{m.h0_km_s_mpc, m.omega_m, m.omega_r, m.omega_b,
                        m.omega_gamma, r.z_drag};
  for (double x : values)
    if (!std::isfinite(x))
      return numerics::Status::nonfinite_input;
  if (!(m.h0_km_s_mpc > 0) || m.omega_m < 0 || !(m.omega_r > 0) ||
      m.omega_b < 0 || m.omega_b > m.omega_m || !(m.omega_gamma > 0) ||
      m.omega_gamma > m.omega_r || r.z_drag < 0 || r.drag_origin.empty())
    return numerics::Status::outside_domain;
  // Subtract the larger fraction first. When it exceeds 1/2, 1-largest
  // is exact for binary64 inputs (Sterbenz); otherwise closure is safely
  // below 1. This rejects 1+minpositive in either order even when a sum rounds1.
  const long double largest = std::max(m.omega_m, m.omega_r);
  const long double smaller = std::min(m.omega_m, m.omega_r);
  const long double remaining = 1.L - largest;
  if (largest > 1 || smaller > remaining)
    return numerics::Status::outside_domain;
  s = {(long double)m.omega_r, (long double)m.omega_m, remaining - smaller,
       (3.L * m.omega_b) / (4.L * m.omega_gamma),
       1.L / (1.L + (long double)r.z_drag)};
  return numerics::Status::ok;
}
double integrand(double t, const void *context) {
  const auto &s = *static_cast<const State *>(context);
  const long double a = s.a_drag * (long double)t;
  const long double a2 = a * a;
  const long double denominator =
      std::sqrt(s.radiation + s.matter * a + s.lambda * a2 * a2) *
      std::sqrt(1.L + s.baryon_loading * a);
  return static_cast<double>(1.L / denominator);
}
} // namespace

std::optional<std::size_t> sound_horizon_payload_bound(
    std::size_t count, std::size_t origins) noexcept {
  if (!count) return origins == 0 ? std::optional<std::size_t>{0} : std::nullopt;
  irred::detail::PayloadAccounting bytes(sizeof(SoundHorizonBatch));
  // Rows reserve once; temporary construction and copied origins coexist.
  bytes.add(count, 2 * sizeof(SoundHorizonRow));
  bytes.add(origins, 2);
  return bytes.result();
}
std::optional<std::size_t> sound_horizon_payload_bound(
    std::span<const SoundHorizonRequest> input) noexcept {
  std::size_t origins = 0;
  for (const auto &r : input) {
    if (!irred::detail::checked_payload_add(origins, r.drag_origin.size(), 1) ||
        !irred::detail::checked_payload_add(origins, 1, 1)) return {};
  }
  return sound_horizon_payload_bound(input.size(), origins);
}

SoundHorizonBatch evaluate_sound_horizon(
    std::span<const SoundHorizonRequest> input, SoundHorizonPolicy p) {
  SoundHorizonBatch out;
  // Qualified arithmetic uses binary64 callbacks and at least 64-bit wide
  // mantissas/exponents for scaled coefficients and physical projections.
  if (std::numeric_limits<long double>::digits < 64 ||
      std::numeric_limits<long double>::max_exponent < 16384 ||
      std::fegetround() != FE_TONEAREST) return out;
  if (!std::isfinite(p.absolute_tolerance_mpc) ||
      !std::isfinite(p.relative_tolerance) || p.absolute_tolerance_mpc < 0 ||
      p.relative_tolerance < 0 ||
      (p.absolute_tolerance_mpc == 0 && p.relative_tolerance == 0) ||
      p.maximum_depth > 60 || input.size() > 65536)
    return out;
  if (input.empty()) {
    out.status = numerics::Status::ok;
    return out;
  }
  const auto bytes = sound_horizon_payload_bound(input);
  if (!bytes)
    return out;
  if (input.size() > p.maximum_points || *bytes > p.maximum_native_bytes) {
    out.status = numerics::Status::work_limit;
    return out;
  }
  out.rows.reserve(input.size());
  out.status = numerics::Status::ok;
  for (const auto &r : input) {
    SoundHorizonRow attempted;
    attempted.source = r;
    out.rows.push_back(std::move(attempted));
    auto &row = out.rows.back();
    State s{};
    row.status = state(r, s);
    if (row.status != numerics::Status::ok)
      continue;
    const std::size_t remaining = p.maximum_total_callbacks - out.callbacks;
    const std::size_t available = std::min(p.maximum_callbacks_per_point, remaining);
    if (available < 3) {
      row.status = numerics::Status::work_limit;
      continue;
    }
    const auto scale = detail::prepare_flat_scale(r.model.h0_km_s_mpc);
    const long double prefactor = scale.distance_mpc * s.a_drag / std::sqrt(3.L);
    const long double converted = (long double)p.absolute_tolerance_mpc / prefactor;
    const double absolute = static_cast<double>(std::min(
        converted, (long double)std::numeric_limits<double>::max()));
    if (absolute == 0 && p.relative_tolerance == 0) {
      row.status = numerics::Status::conditioning_budget_exceeded;
      continue;
    }
    const auto q = numerics::integrate(
        integrand, &s, 0, 1,
        {absolute, p.relative_tolerance, available, p.maximum_depth});
    row.callbacks = q.evaluations;
    out.callbacks += q.evaluations;
    row.status = q.status;
    if (q.status != numerics::Status::ok)
      continue;
    const long double value = prefactor * (long double)q.value;
    const long double estimate = prefactor * (long double)q.error_estimate;
    if (!detail::physical_representable(value) || !(value > 0) ||
        !detail::representable(estimate)) {
      row.status = numerics::Status::outside_domain;
      continue;
    }
    const double rounded = static_cast<double>(value);
    const long double integral_rounding =
        ((long double)std::nextafter(q.value, std::numeric_limits<double>::infinity()) -
         (long double)q.value) / 2;
    // Positive callback values have no cancellation in the integral. Reserve
    // a conservative arithmetic admission floor: 32 binary64 eps for callback
    // conversion/weighted integral storage, and 128 wide eps for coefficient,
    // square-root and scale arithmetic. This is an admission allowance, not
    // a rigorous libm or quadrature bound; independent refinements remain gates.
    const long double arithmetic_floor = std::abs(value) *
        (32.L * std::numeric_limits<double>::epsilon() +
         128.L * std::numeric_limits<long double>::epsilon());
    const long double diagnostic = estimate + arithmetic_floor + prefactor * integral_rounding +
                                   std::abs(value - (long double)rounded);
    const long double budget = (long double)p.absolute_tolerance_mpc +
                               (long double)p.relative_tolerance * std::abs(value);
    if (!detail::representable(diagnostic) || !(diagnostic <= budget)) {
      row.status = numerics::Status::conditioning_budget_exceeded;
      continue;
    }
    row.sound_horizon_mpc = rounded;
    row.error_estimate_mpc = static_cast<double>(diagnostic);
  }
  return out;
}
} // namespace irred::cosmology
