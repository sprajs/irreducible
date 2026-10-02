#include "irred/bao_thermal.hpp"
#include "bao_density_projection.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace irred::bao {
std::optional<size_t>
thermal_density_payload_bound(size_t m, size_t n, size_t strings,
                              size_t species, size_t maximum_species,
                              unsigned requested) noexcept {
  if (!requested || (requested & ~7u) || (n && m > SIZE_MAX / n))
    return {};
  irred::detail::PayloadAccounting b(sizeof(ThermalDensityBatch));
  b.add(m, sizeof(ThermalDensitySlot));
  b.add(strings, 1);
  b.add(species, sizeof(cosmology::ThermalPhysicalSpecies));
  b.add(m * n, (bool(requested & 2u) + bool(requested & 4u)) * sizeof(double));
  b.add(n, 12 * sizeof(double) + 16 * sizeof(long double));
  const auto provider =
      cosmology::thermal_observables_payload_bound(n, maximum_species, strings);
  if (!provider)
    return {};
  b.add(*provider, 1);
  return b.result();
}

ThermalDensityBatch PreparedDensity::evaluate_thermal(
    std::span<const cosmology::ThermalObservableRequest> points,
    ThermalDensityPolicy p) const {
  ThermalDensityBatch out;
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
  irred::detail::PayloadAccounting origins(0), copied_origins(0), species(0);
  size_t maximum_species = 0;
  for (const auto &x : points) {
    for (const auto *origin : {&x.drag_origin, &x.source_origin}) {
      origins.add(origin->size(), 1);
      if (origin->size() == SIZE_MAX)
        copied_origins.add(SIZE_MAX, 2);
      else
        copied_origins.add(std::max(size_t(32), origin->size() + 1), 1);
    }
    species.add(x.model.species.size(), 1);
    maximum_species = std::max(maximum_species, x.model.species.size());
  }
  if (maximum_species > 16) {
    out.numerical_status = numerics::Status::outside_domain;
    return out;
  }
  const auto strings = origins.result();
  const auto copied = copied_origins.result();
  const auto count = species.result();
  const auto bytes =
      copied && count
          ? thermal_density_payload_bound(points.size(), n, *copied, *count,
                                          maximum_species, p.requested)
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
    local.thermal.maximum_total_callbacks =
        std::min(local.thermal.maximum_total_callbacks, remaining);
    // Mapping and normalization occur once for this model, before its coarse
    // ratio evaluation; normalization failures still charge actual work.
    const auto owner = cosmology::prepare_thermal_observables(point, local);
    s.preparation_status = owner.status();
    s.preparation_callbacks = owner.background().preparation_callbacks();
    s.momentum_callbacks = s.preparation_callbacks;
    s.callbacks = s.preparation_callbacks;
    out.preparation_callbacks += s.preparation_callbacks;
    out.momentum_callbacks += s.preparation_callbacks;
    out.callbacks += s.preparation_callbacks;
    remaining -= s.preparation_callbacks;
    auto cause = owner.status();
    cosmology::ThermalObservableBatch prediction;
    if (cause == numerics::Status::ok) {
      local.maximum_total_callbacks = std::min(
          local.maximum_total_callbacks - s.preparation_callbacks, remaining);
      local.thermal.maximum_total_callbacks -= s.preparation_callbacks;
      prediction = owner.evaluate(z, mask, local);
      s.outer_callbacks = prediction.outer_callbacks;
      s.momentum_callbacks += prediction.momentum_callbacks;
      s.callbacks += prediction.callbacks;
      out.outer_callbacks += prediction.outer_callbacks;
      out.momentum_callbacks += prediction.momentum_callbacks;
      out.callbacks += prediction.callbacks;
      remaining -= prediction.callbacks;
      cause = prediction.status;
    }
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
