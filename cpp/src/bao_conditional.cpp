#include "irred/bao_conditional.hpp"
#include "bao_density_projection.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace irred::bao {
std::optional<size_t>
conditional_density_payload_bound(size_t m, size_t n, size_t strings,
                                  unsigned requested) noexcept {
  if (!requested || (requested & ~7u) || (n && m > SIZE_MAX / n))
    return {};
  irred::detail::PayloadAccounting b(sizeof(ConditionalDensityBatch));
  b.add(m, sizeof(ConditionalDensitySlot));
  b.add(strings, 1);
  b.add(m, 32); // copied-origin SSO allowance
  b.add(m * n, (bool(requested & 2u) + bool(requested & 4u)) * sizeof(double));
  b.add(n, 12 * sizeof(double) +
               16 * sizeof(long double)); // serial prediction/solve
  // At most the sum of all origins is used by any serial provider invocation.
  const auto provider = cosmology::early_late_payload_bound(n, strings);
  if (!provider)
    return {};
  b.add(*provider, 1);
  return b.result();
}

ConditionalDensityBatch PreparedDensity::evaluate_conditional(
    std::span<const cosmology::SoundHorizonRequest> points,
    ConditionalDensityPolicy p) const {
  ConditionalDensityBatch out;
  if (status_ != statistics::DensityStatus::finite) {
    out.status = status_;
    out.numerical_status = numerical_status_;
    return out;
  }
  if (!p.requested || (p.requested & ~7u) ||
      !std::isfinite(p.maximum_forward_sensitivity) ||
      p.maximum_forward_sensitivity <= 0 ||
      !std::isfinite(p.maximum_projection_log_density_error) ||
      p.maximum_projection_log_density_error <= 0)
    return out;
  if (p.arithmetic != gaussian_.arithmetic()) {
    out.status = statistics::DensityStatus::incompatible_metadata;
    return out;
  }
  const size_t n = source_.queries.size();
  if (points.size() > p.maximum_models || n > p.maximum_queries) {
    out.numerical_status = numerics::Status::work_limit;
    return out;
  }
  irred::detail::PayloadAccounting origins(0), copied_origins(0);
  for (const auto &x : points) {
    origins.add(x.drag_origin.size(), 1);
    if (x.drag_origin.size() == SIZE_MAX)
      copied_origins.add(SIZE_MAX, 2);
    else
      copied_origins.add(std::max(size_t(32), x.drag_origin.size() + 1), 1);
  }
  const auto strings = origins.result();
  const auto copied = copied_origins.result();
  const auto bytes = copied ? conditional_density_payload_bound(
                                  points.size(), n, *copied, p.requested)
                            : std::nullopt;
  if (!strings || *strings > p.maximum_string_bytes || !bytes ||
      *bytes > p.maximum_native_bytes) {
    out.numerical_status = numerics::Status::work_limit;
    return out;
  }
  std::vector<double> z;
  z.reserve(n);
  unsigned mask = 0;
  for (const auto &q : source_.queries) {
    z.push_back(q.z);
    mask |= cosmology::early_late_mask(detail::ratio_output(q.observable));
  }
  out.slots.reserve(points.size());
  out.status = statistics::DensityStatus::finite;
  out.numerical_status = numerics::Status::ok;
  long double inverse_norm_estimate = 0;
  if (p.requested & 1u) {
    inverse_norm_estimate = detail::inverse_covariance_norm_estimate(
        gaussian_.factor_, source_.covariance, n);
  }
  size_t remaining = p.maximum_total_callbacks;
  for (const auto &point : points) {
    out.slots.emplace_back();
    auto &s = out.slots.back();
    s.source = point;
    auto local = p.predictions;
    local.maximum_total_callbacks =
        std::min(local.maximum_total_callbacks, remaining);
    local.sound.maximum_total_callbacks =
        std::min(local.sound.maximum_total_callbacks, remaining);
    const auto prediction =
        cosmology::evaluate_early_late(point, z, mask, local);
    s.callbacks = prediction.callbacks;
    out.callbacks += s.callbacks;
    remaining -= s.callbacks;
    auto cause = prediction.status;
    std::vector<double> mu, r;
    std::vector<long double> eps;
    mu.reserve(n);
    if (p.requested & 5u)
      r.reserve(n);
    if (p.requested & 1u)
      eps.reserve(n);
    if (cause == numerics::Status::ok) {
      for (size_t j = 0; j < n; ++j) {
        const auto &v = prediction.rows[j].outputs[static_cast<unsigned>(
            detail::ratio_output(source_.queries[j].observable))];
        if (v.status != numerics::Status::ok || !v.value) {
          cause = v.status;
          break;
        }
        mu.push_back(*v.value);
        if (p.requested & 1u)
          eps.push_back(v.error_estimate);
      }
    }
    if (p.requested & 2u)
      detail::output_state(s.predictions_state, cause);
    if (cause == numerics::Status::ok && (p.requested & 5u)) {
      cause = detail::subtract_observations(source_.observed, mu, r, eps);
    }
    if (p.requested & 4u)
      detail::output_state(s.residuals_state, cause);
    if (cause == numerics::Status::ok && (p.requested & 1u)) {
      cause = detail::project_density(
          gaussian_.factor_, gaussian_, source_.ordered_ids, r, eps,
          inverse_norm_estimate, p.maximum_forward_sensitivity,
          p.maximum_projection_log_density_error,
          s.projection_log_density_error_estimate, s.result);
    }
    if (p.requested & 1u)
      detail::output_state(s.density_state, cause);
    s.numerical_status = cause;
    if ((p.requested & 2u) && mu.size() == n)
      s.predictions = std::move(mu);
    if ((p.requested & 4u) && r.size() == n)
      s.residuals = std::move(r);
  }
  return out;
}
} // namespace irred::bao
