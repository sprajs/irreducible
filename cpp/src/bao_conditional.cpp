#include "irred/bao_conditional.hpp"
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
  detail::PayloadAccounting b(sizeof(ConditionalDensityBatch));
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
namespace {
void state(OutputState &s, numerics::Status cause) {
  s.availability = cause == numerics::Status::ok
                       ? cosmology::Availability::available
                       : cosmology::Availability::failed;
  s.numerical_status = cause;
  s.status = cause == numerics::Status::ok ? cosmology::Status::ok
             : cause == numerics::Status::work_limit
                 ? cosmology::Status::work_limit
                 : cosmology::Status::numerical_failure;
}
cosmology::EarlyLateOutput output(Observable x) {
  switch (x) {
  case Observable::transverse_over_ruler:
    return cosmology::EarlyLateOutput::dm_over_rs;
  case Observable::hubble_over_ruler:
    return cosmology::EarlyLateOutput::dh_over_rs;
  case Observable::volume_over_ruler:
    return cosmology::EarlyLateOutput::dv_over_rs;
  }
  return cosmology::EarlyLateOutput::count;
}
} // namespace
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
  detail::PayloadAccounting origins(0), copied_origins(0);
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
    mask |= cosmology::early_late_mask(output(q.observable));
  }
  out.slots.reserve(points.size());
  out.status = statistics::DensityStatus::finite;
  out.numerical_status = numerics::Status::ok;
  long double inverse_norm_estimate = 0;
  if (p.requested & 1u) {
    long double covariance_norm = 0;
    for (size_t j = 0; j < n; ++j) {
      long double row = 0;
      for (size_t k = 0; k < n; ++k)
        row += std::abs((long double)source_.covariance[j * n + k]);
      covariance_norm = std::max(covariance_norm, row);
    }
    inverse_norm_estimate =
        gaussian_.factor_.condition_estimate_inf() / covariance_norm;
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
            output(source_.queries[j].observable))];
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
      state(s.predictions_state, cause);
    if (cause == numerics::Status::ok && (p.requested & 5u)) {
      for (size_t j = 0; j < n; ++j) {
        const long double wide = (long double)source_.observed[j] - mu[j];
        const double cast = static_cast<double>(wide);
        if (!std::isfinite(cast) ||
            (cast != 0 && std::fpclassify(cast) != FP_NORMAL) ||
            (cast == 0 && wide != 0)) {
          cause = numerics::Status::overflow;
          break;
        }
        r.push_back(cast);
        if (p.requested & 1u)
          eps[j] += std::abs(wide - cast) +
                    std::numeric_limits<long double>::epsilon() *
                        (std::abs((long double)source_.observed[j]) +
                         std::abs((long double)mu[j]));
      }
    }
    if (p.requested & 4u)
      state(s.residuals_state, cause);
    if (cause == numerics::Status::ok && (p.requested & 1u)) {
      const auto solved =
          numerics::solve(gaussian_.factor_, r, p.maximum_forward_sensitivity);
      cause = solved.status;
      if (cause == numerics::Status::ok) {
        long double ainf = 0, einf = 0, e1 = 0;
        for (size_t j = 0; j < n; ++j) {
          ainf = std::max(ainf, std::abs((long double)solved.value[j]));
          einf = std::max(einf, eps[j]);
          e1 += eps[j];
        }
        const long double estimate =
            ainf * e1 + inverse_norm_estimate * einf * e1 / 2;
        const double reported = static_cast<double>(estimate);
        if (!std::isfinite(estimate) || estimate < 0 ||
            !std::isfinite(reported))
          cause = numerics::Status::overflow;
        else if (estimate != 0 &&
                 (reported == 0 || std::fpclassify(reported) != FP_NORMAL))
          cause = numerics::Status::outside_domain;
        else {
          s.projection_log_density_error_estimate = reported;
          if (estimate > p.maximum_projection_log_density_error)
            cause = numerics::Status::conditioning_budget_exceeded;
        }
      }
      if (cause == numerics::Status::ok) {
        s.result = gaussian_.evaluate(r, source_.ordered_ids,
                                      p.maximum_forward_sensitivity);
        cause = s.result->density.numerical_status;
      }
    }
    if (p.requested & 1u)
      state(s.density_state, cause);
    s.numerical_status = cause;
    if ((p.requested & 2u) && mu.size() == n)
      s.predictions = std::move(mu);
    if ((p.requested & 4u) && r.size() == n)
      s.residuals = std::move(r);
  }
  return out;
}
} // namespace irred::bao
