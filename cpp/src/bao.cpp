#include "irred/bao.hpp"
#include "payload_accounting.hpp"
#include "piecewise_radial.hpp"
#include <cmath>
#include <limits>
#include <utility>
namespace irred::bao {
PreparedDensity::PreparedDensity(PreparedDensity &&other) noexcept
    : source_(std::move(other.source_)), gaussian_(std::move(other.gaussian_)),
      status_(std::exchange(other.status_,
                            statistics::DensityStatus::invalid_input)),
      numerical_status_(std::exchange(other.numerical_status_,
                                      numerics::Status::invalid_input)) {}
PreparedDensity &PreparedDensity::operator=(PreparedDensity &&other) noexcept {
  if (this != &other) {
    source_ = std::move(other.source_);
    gaussian_ = std::move(other.gaussian_);
    status_ =
        std::exchange(other.status_, statistics::DensityStatus::invalid_input);
    numerical_status_ =
        std::exchange(other.numerical_status_, numerics::Status::invalid_input);
  }
  return *this;
}

std::optional<size_t> PreparedDensity::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  const auto source = retained_source_payload_bound(source_);
  if (!source)
    return {};
  b.add(*source, 1);
  b.embedded(gaussian_.retained_payload_bound(), sizeof(gaussian_));
  return b.result();
}
namespace {
bool normal_or_zero(long double x) {
  return std::isfinite(x) &&
         std::abs(x) <= std::numeric_limits<double>::max() &&
         (x == 0 ||
          ((double)x != 0 && std::fpclassify((double)x) == FP_NORMAL));
}
bool known(Observable o) {
  return o == Observable::transverse_over_ruler ||
         o == Observable::hubble_over_ruler ||
         o == Observable::volume_over_ruler;
}
template <class P> bool valid_policy(const P &p) {
  return std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0 &&
         (p.arithmetic == numerics::Arithmetic::binary64_legacy_v1 ||
          p.arithmetic == numerics::Arithmetic::longdouble_cpu_v1);
}
bool fits(size_t count, size_t width, size_t &remaining) {
  if (count > remaining / width)
    return false;
  remaining -= count * width;
  return true;
}
long double project(long double scale, long double E, long double I, Query q) {
  const auto dm = scale * I;
  const auto dh = scale / E;
  switch (q.observable) {
  case Observable::transverse_over_ruler:
    return dm;
  case Observable::hubble_over_ruler:
    return dh;
  case Observable::volume_over_ruler:
    return std::cbrt((long double)q.z * dm * dm * dh);
  }
  return 0;
}
} // namespace
Batch evaluate(const cosmology::Expansion &background, Ruler ruler,
               std::span<const Query> queries, Policy policy) {
  Batch out;
  const auto background_payload =
      cosmology::Expansion::workspace_payload_bound(queries.size());
  auto available_bytes = policy.maximum_native_bytes;
  if (queries.size() > policy.maximum_queries || !background_payload ||
      !fits(*background_payload, 1, available_bytes) ||
      !fits(queries.size(),
            sizeof(Slot) + sizeof(cosmology::Request) + sizeof(size_t),
            available_bytes)) {
    out.status = cosmology::Status::work_limit;
    return out;
  }
  if (!std::isfinite(ruler.h0_rd_km_s))
    return out;
  if (ruler.h0_rd_km_s < 5000 || ruler.h0_rd_km_s > 15000) {
    out.status = cosmology::Status::unsupported_domain;
    return out;
  }
  if (background.status() != cosmology::Status::ok) {
    out.status = background.status();
    return out;
  }
  std::vector<cosmology::Request> requests;
  std::vector<size_t> indices;
  requests.reserve(queries.size());
  indices.reserve(queries.size());
  out.slots.resize(queries.size());
  for (size_t i = 0; i < queries.size(); ++i) {
    out.slots[i].source = queries[i];
    if (!known(queries[i].observable))
      continue;
    requests.emplace_back(
        queries[i].z,
        (uint32_t)(queries[i].observable == Observable::hubble_over_ruler
                       ? cosmology::Observable::expansion
                       : cosmology::Observable::radial));
    indices.push_back(i);
  }
  const auto computed = background.evaluate(requests, policy.background);
  out.status = computed.status;
  out.work = computed.work;
  if (computed.status != cosmology::Status::ok) {
    out.slots.clear();
    return out;
  }
  const auto scale =
      (long double)irred::speed_of_light_m_per_s / 1000 / ruler.h0_rd_km_s;
  for (size_t j = 0; j < indices.size(); ++j) {
    auto &s = out.slots[indices[j]];
    const auto &r = computed.slots[j];
    s.node_index = r.node_index;
    long double E = 0, I = 0;
    if (s.source.observable == Observable::hubble_over_ruler) {
      s.status = r.expansion.status;
      s.numerical_status = r.expansion.numerical_status;
      if (!r.expansion.value)
        continue;
      E = cosmology::detail::RadialAccess::expansion(*r.expansion.value);
    } else {
      s.status = r.radial.status;
      s.numerical_status = r.radial.numerical_status;
      if (!r.radial.value)
        continue;
      E = cosmology::detail::RadialAccess::expansion(*r.radial.value);
      I = cosmology::detail::RadialAccess::integral(*r.radial.value);
    }
    const auto v = project(scale, E, I, s.source);
    if (!normal_or_zero(v)) {
      s.status = cosmology::Status::numerical_failure;
      s.numerical_status = std::abs(v) > std::numeric_limits<double>::max()
                               ? numerics::Status::overflow
                               : numerics::Status::outside_domain;
      continue;
    }
    s.value = (double)v;
    s.status = cosmology::Status::ok;
    s.numerical_status = numerics::Status::ok;
  }
  return out;
}
std::optional<size_t>
retained_source_payload_bound(const DensityInput &input) noexcept {
  size_t remaining = std::numeric_limits<size_t>::max();
  if (!fits(input.queries.capacity(), sizeof(Query), remaining) ||
      !fits(input.observed.capacity(), sizeof(double), remaining) ||
      !fits(input.covariance.capacity(), sizeof(double), remaining) ||
      !fits(input.ordered_ids.capacity(), sizeof(std::string), remaining))
    return {};
  auto string = [&](const std::string &value) {
    return value.capacity() != std::numeric_limits<size_t>::max() &&
           fits(value.capacity() + 1, 1, remaining);
  };
  for (const auto &id : input.ordered_ids)
    if (!string(id))
      return {};
  for (const auto *value :
       {&input.table_identity, &input.covariance_identity,
        &input.ordering_provenance, &input.calibration_provenance,
        &input.dependence_provenance})
    if (!string(*value))
      return {};
  return std::numeric_limits<size_t>::max() - remaining;
}
std::optional<size_t> preparation_payload_bound(
    size_t n, size_t source_payload, size_t identity_payload) noexcept {
  if (n && n > SIZE_MAX / n)
    return {};
  detail::PayloadAccounting b(sizeof(PreparedDensity));
  b.add(source_payload, 1);
  b.add(identity_payload, 1);
  b.add(n, sizeof(std::string));
  b.add(n * n, 2 * sizeof(double) + sizeof(long double));
  b.add(n, 8 * sizeof(long double) + 4 * sizeof(double));
  // Conservatively admit compiled metadata strings, including SSO storage.
  b.add(1, 4096);
  return b.result();
}
PreparedDensity prepare_density(DensityInput input, PreparationPolicy policy) {
  PreparedDensity out;
  const auto n = input.queries.size();
  if (!valid_policy(policy))
    return out;
  if (n > policy.maximum_queries || n == 0 ||
      n > policy.maximum_matrix_elements / n) {
    out.numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  if (input.observed.size() != n || input.ordered_ids.size() != n ||
      input.covariance.size() != n * n)
    return out;
  const auto source_payload = retained_source_payload_bound(input);
  detail::PayloadAccounting copied(0);
  for (const auto &id : input.ordered_ids)
    copied.string(id);
  for (const auto *text :
       {&input.table_identity, &input.covariance_identity,
        &input.ordering_provenance, &input.calibration_provenance,
        &input.dependence_provenance})
    copied.string(*text);
  const auto identity_payload = copied.result();
  const auto peak = source_payload && identity_payload
                        ? preparation_payload_bound(n, *source_payload,
                                                    *identity_payload)
                        : std::nullopt;
  if (!peak || *peak > policy.maximum_native_bytes) {
    out.numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  auto string_bytes = policy.maximum_string_bytes;
  for (const auto *text :
       {&input.table_identity, &input.covariance_identity,
        &input.ordering_provenance, &input.calibration_provenance,
        &input.dependence_provenance}) {
    if (text->empty())
      return out;
    if (!fits(text->size(), 1, string_bytes)) {
      out.numerical_status_ = numerics::Status::work_limit;
      return out;
    }
  }
  if ((input.role != RowRole::released_fitted_distance_summary &&
       input.role != RowRole::synthetic_control) ||
      input.covariance_unit != CovarianceUnit::dimensionless_ratio_squared)
    return out;
  for (std::size_t i = 0; i < n; ++i) {
    if (input.ordered_ids[i].empty() || !known(input.queries[i].observable) ||
        !std::isfinite(input.queries[i].z) || input.queries[i].z < 0 ||
        input.queries[i].z > 5 || !std::isfinite(input.observed[i]))
      return out;
    if (!fits(input.ordered_ids[i].size(), 1, string_bytes)) {
      out.numerical_status_ = numerics::Status::work_limit;
      return out;
    }
    for (std::size_t j = 0; j < i; ++j)
      if (input.ordered_ids[i] == input.ordered_ids[j])
        return out;
  }
  statistics::Metadata metadata;
  metadata.ordered_ids = input.ordered_ids;
  metadata.measure =
      "product of dimensionless BAO released distance-ratio coordinates";
  metadata.table_identity = input.table_identity;
  metadata.uncertainty_identity = input.covariance_identity;
  metadata.ordering_provenance = input.ordering_provenance;
  metadata.calibration_provenance = input.calibration_provenance;
  metadata.dependence_provenance = input.dependence_provenance;
  metadata.source_semantics =
      input.role == RowRole::synthetic_control
          ? "synthetic BAO control; covariance dimensionless ratio squared"
          : "released fitted BAO distance summaries; covariance dimensionless "
            "ratio squared";
  metadata.input_matrix_convention =
      "full source covariance in declared exact row order";
  out.gaussian_ = statistics::prepare_gaussian(
      input.covariance, statistics::MatrixKind::covariance, std::move(metadata),
      policy.maximum_matrix_elements, policy.maximum_forward_sensitivity,
      policy.arithmetic);
  out.status_ = out.gaussian_.status();
  out.numerical_status_ = out.gaussian_.numerical_status();
  out.source_ = std::move(input);
  return out;
}
DensityBatch PreparedDensity::evaluate(std::span<const ModelPoint> models,
                                       DensityPolicy p) const {
  DensityBatch out;
  if (status_ != statistics::DensityStatus::finite) {
    out.status = status_;
    out.numerical_status = numerical_status_;
    return out;
  }
  if (!valid_policy(p) || !p.requested || (p.requested & ~7u))
    return out;
  if (p.arithmetic != gaussian_.arithmetic()) {
    out.status = statistics::DensityStatus::incompatible_metadata;
    return out;
  }
  const auto n = source_.queries.size();
  auto bytes = p.maximum_native_bytes;
  const auto background_payload =
      cosmology::Expansion::workspace_payload_bound(n);
  if (models.size() > p.maximum_models ||
      (n && models.size() > std::numeric_limits<size_t>::max() / n) ||
      !background_payload || !fits(models.size(), sizeof(DensitySlot), bytes) ||
      !fits(models.size() * n,
            (bool(p.requested & 2u) + bool(p.requested & 4u)) * sizeof(double) +
                sizeof(size_t),
            bytes) ||
      !fits(*background_payload, 1, bytes) ||
      !fits(n,
            sizeof(Slot) + sizeof(cosmology::Request) + sizeof(size_t) +
                8 * sizeof(double) + 8 * sizeof(long double),
            bytes)) {
    out.numerical_status = numerics::Status::work_limit;
    return out;
  }
  out.slots.reserve(models.size());
  out.status = statistics::DensityStatus::finite;
  out.numerical_status = numerics::Status::ok;
  auto callbacks = p.observables.background.maximum_callbacks,
       segments = p.observables.background.maximum_segment_visits;
  for (const auto &point : models) {
    out.slots.emplace_back(point);
    auto &s = out.slots.back();
    auto background = cosmology::prepare(point.expansion, point.geometry);
    s.model_id = background.model_id();
    auto local = p.observables;
    local.background.maximum_callbacks = callbacks;
    local.background.maximum_segment_visits = segments;
    const auto predicted =
        bao::evaluate(background, point.ruler, source_.queries, local);
    s.background_status = predicted.status;
    s.work = predicted.work;
    callbacks -= predicted.work.callbacks;
    segments -= predicted.work.segment_visits;
    out.work.callbacks += predicted.work.callbacks;
    out.work.segment_visits += predicted.work.segment_visits;
    auto set = [](OutputState &state, bool ok, numerics::Status cause) {
      state.availability = ok ? cosmology::Availability::available
                              : cosmology::Availability::failed;
      state.status = ok ? cosmology::Status::ok
                     : cause == numerics::Status::work_limit
                         ? cosmology::Status::work_limit
                         : cosmology::Status::numerical_failure;
      state.numerical_status = cause;
    };
    auto fail_density = [&](numerics::Status cause) {
      if (p.requested & 1u) {
        set(s.density_state, false, cause);
        s.result.emplace();
        s.result->density.status = statistics::DensityStatus::numerical_failure;
        s.result->density.numerical_status = cause;
      }
    };
    if (predicted.status != cosmology::Status::ok) {
      s.numerical_status =
          predicted.status == cosmology::Status::work_limit
              ? numerics::Status::work_limit
          : predicted.status == cosmology::Status::unsupported_domain
              ? numerics::Status::outside_domain
              : numerics::Status::invalid_input;
      if (p.requested & 2u)
        set(s.predictions_state, false, s.numerical_status);
      if (p.requested & 4u)
        set(s.residuals_state, false, s.numerical_status);
      fail_density(s.numerical_status);
      continue;
    }
    std::vector<double> predictions, residuals;
    predictions.reserve(n);
    s.node_indices.reserve(n);
    bool finite = true;
    for (size_t j = 0; j < n; ++j) {
      const auto &r = predicted.slots[j];
      if (!r.value) {
        s.background_status = r.status;
        s.numerical_status = r.numerical_status;
        finite = false;
        break;
      }
      predictions.push_back(*r.value);
      s.node_indices.push_back(*r.node_index);
    }
    if (p.requested & 2u)
      set(s.predictions_state, finite,
          finite ? numerics::Status::ok : s.numerical_status);
    if (!finite) {
      s.node_indices.clear();
      if (p.requested & 4u)
        set(s.residuals_state, false, s.numerical_status);
      fail_density(s.numerical_status);
      continue;
    }
    s.has_predictions = true;
    s.numerical_status = numerics::Status::ok;
    if (p.requested & 5u) {
      residuals.reserve(n);
      for (size_t j = 0; j < n; ++j) {
        const auto residual = (long double)source_.observed[j] - predictions[j];
        if (!normal_or_zero(residual)) {
          s.numerical_status =
              std::abs(residual) > std::numeric_limits<double>::max()
                  ? numerics::Status::overflow
                  : numerics::Status::outside_domain;
          finite = false;
          break;
        }
        residuals.push_back((double)residual);
      }
      if (p.requested & 4u)
        set(s.residuals_state, finite,
            finite ? numerics::Status::ok : s.numerical_status);
      if (!finite) {
        if (p.requested & 2u)
          s.predictions = std::move(predictions);
        fail_density(s.numerical_status);
        continue;
      }
    }
    if (p.requested & 1u) {
      s.result = gaussian_.evaluate(residuals, source_.ordered_ids,
                                    p.maximum_forward_sensitivity);
      s.numerical_status = s.result->density.numerical_status;
      set(s.density_state,
          s.result->density.status == statistics::DensityStatus::finite,
          s.numerical_status);
    }
    if (p.requested & 2u)
      s.predictions = std::move(predictions);
    if (p.requested & 4u)
      s.residuals = std::move(residuals);
  }
  return out;
}
} // namespace irred::bao
