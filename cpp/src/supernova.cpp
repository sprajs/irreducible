#include "irred/supernova.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace irred::supernova {
std::optional<size_t> Consumer::retained_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(*this));
  b.vector(selected_.source_indices);
  b.strings(selected_.ordered_ids);
  b.vector(selected_.coordinates);
  b.vector(observed_);
  b.embedded(profile_.retained_payload_bound(), sizeof(profile_));
  return b.result();
}
namespace {
bool normal(long double v) {
  return std::isfinite(v) &&
         std::abs(v) <= std::numeric_limits<double>::max() &&
         (v == 0 ||
          ((double)v != 0 && std::fpclassify((double)v) == FP_NORMAL));
}
constexpr auto score_bit = (std::uint32_t)Output::score,
               geometry_bit = (std::uint32_t)Output::geometric_shape,
               effect_bit = (std::uint32_t)Output::magnitude_effect,
               corrected_bit = (std::uint32_t)Output::corrected_residuals,
               profiled_bit = (std::uint32_t)Output::profiled_residuals,
               diagnostic_bit = (std::uint32_t)Output::diagnostics;
} // namespace
std::optional<size_t>
preparation_payload_bound(const observations::Prepared &source, size_t selected,
                          numerics::Arithmetic arithmetic) noexcept {
  if (arithmetic > numerics::Arithmetic::longdouble_cpu_v1 ||
      (selected && selected > SIZE_MAX / selected))
    return {};
  detail::PayloadAccounting b(sizeof(Consumer));
  // Supported vector implementation: push-grown complement metadata capacity
  // below twice raw rows; fixed metadata text envelope 4096 bytes.
  const auto &s = source.source();
  b.add(selected,
        sizeof(size_t) + sizeof(MagnitudeCoordinate) + sizeof(std::string));
  b.add(s.values.size(), 4 * sizeof(std::string) + sizeof(uint8_t));
  for (const auto &id : s.measurement_ids) {
    if (id.capacity() == SIZE_MAX)
      return {};
    b.add(id.capacity() + 1, 3);
  }
  for (const auto *t :
       {&s.table_sha256, &s.uncertainty_sha256, &s.ordering_provenance,
        &s.calibration_provenance, &s.dependence_provenance})
    b.string(*t);
  b.add(1, 4096);
  // Simultaneous local selected block, Gaussian covariance, canonical factor
  // original, lower storage (both policies conservatively admitted), plus
  // serial basis solve/profile buffers. No inverse matrix is retained here.
  b.add(selected * selected, 4 * sizeof(double) + sizeof(long double));
  b.add(selected, 16 * sizeof(double) + 16 * sizeof(long double));
  return b.result();
}
Consumer prepare(SelectedMagnitudeSource input, PreparationPolicy p) {
  Consumer out;
  if (!input.source || input.source->status() != observations::Status::ok ||
      p.arithmetic > numerics::Arithmetic::longdouble_cpu_v1 ||
      !std::isfinite(p.maximum_forward_sensitivity) ||
      !(p.maximum_forward_sensitivity > 0))
    return out;
  const auto n = input.source_indices.size();
  if (!n || n > p.maximum_selected_rows || n > p.maximum_matrix_elements / n) {
    out.status_ = Status::work_limit;
    out.preparation_numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  const auto &s = input.source->source();
  if (s.unit != observations::Unit::magnitude ||
      s.uncertainty != observations::Uncertainty::covariance ||
      s.uncertainty_unit != observations::UncertaintyUnit::magnitude_squared ||
      (s.role != observations::Role::observed_measurement &&
       s.role != observations::Role::released_fitted_summary &&
       s.role != observations::Role::synthetic_control) ||
      input.coordinates.size() != n || input.ordered_ids.size() != n) {
    out.status_ = Status::incompatible_metadata;
    return out;
  }
  for (size_t k = 0; k < n; ++k) {
    const auto i = input.source_indices[k];
    const auto &q = input.coordinates[k];
    if (i >= s.values.size() || (k && i <= input.source_indices[k - 1]) ||
        (!s.source_selection.empty() && !s.source_selection[i]) ||
        s.missing[i] || !std::isfinite(s.values[i]) ||
        input.ordered_ids[k] != s.measurement_ids[i] ||
        !std::isfinite(q.z_expansion) || !std::isfinite(q.observer.redshift) ||
        q.observer.redshift <= -1 ||
        (q.observer.convention != cosmology::Convention::released_zhd_zhel &&
         q.observer.convention !=
             cosmology::Convention::geometric_same_redshift) ||
        (q.observer.convention ==
             cosmology::Convention::geometric_same_redshift &&
         q.observer.redshift != q.z_expansion)) {
      out.status_ = Status::incompatible_metadata;
      return out;
    }
  }
  const auto peak = preparation_payload_bound(*input.source, n, p.arithmetic);
  detail::PayloadAccounting required(peak.value_or(0));
  // Account reserve-heavy transferred descriptor capacities above the
  // deterministic n-element envelope before any factor allocation.
  if (input.source_indices.capacity() > n)
    required.add(input.source_indices.capacity() - n, sizeof(size_t));
  if (input.coordinates.capacity() > n)
    required.add(input.coordinates.capacity() - n, sizeof(MagnitudeCoordinate));
  if (input.ordered_ids.capacity() > n)
    required.add(input.ordered_ids.capacity() - n, sizeof(std::string));
  for (size_t k = 0; k < n; ++k) {
    const auto &actual = input.ordered_ids[k];
    const auto &base = s.measurement_ids[input.source_indices[k]];
    if (actual.capacity() > base.capacity())
      required.add(actual.capacity() - base.capacity(), 1);
  }
  const auto payload = required.result();
  if (!peak || !payload || *payload > p.maximum_native_bytes) {
    out.status_ = Status::work_limit;
    out.preparation_numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  // Retain the admitted selection even if factor/profile preparation fails.
  // Moving its descriptor does not duplicate the immutable source or matrix.
  out.selected_ = std::move(input);
  const auto &selected = out.selected_;
  auto g = statistics::prepare_selected_observations(
      *selected.source, selected.source_indices, p.maximum_matrix_elements,
      p.maximum_forward_sensitivity, p.arithmetic);
  out.preparation_status_ = g.status();
  out.preparation_numerical_status_ = g.numerical_status();
  if (g.status() != statistics::DensityStatus::finite) {
    out.status_ = Status::numerical_failure;
    return out;
  }
  out.observed_.reserve(n);
  for (auto i : selected.source_indices)
    out.observed_.push_back(s.values[i]);
  std::vector<double> ones(n, 1);
  out.profile_ = std::move(g).prepare_offset_profile(
      ones, selected.ordered_ids, p.maximum_forward_sensitivity);
  out.preparation_status_ = out.profile_.status();
  out.preparation_numerical_status_ = out.profile_.numerical_status();
  if (out.preparation_status_ != statistics::DensityStatus::finite) {
    out.status_ = Status::numerical_failure;
    return out;
  }
  out.status_ = Status::ok;
  return out;
}
BatchResult Consumer::evaluate_batch(std::span<const ModelPoint> points,
                                     Policy p) const {
  BatchResult out;
  if (!p.requested || (p.requested & ~63u) ||
      p.arithmetic > numerics::Arithmetic::longdouble_cpu_v1 ||
      !std::isfinite(p.maximum_forward_sensitivity) ||
      !(p.maximum_forward_sensitivity > 0))
    return out;
  if (status_ != Status::ok) {
    out.status = status_;
    return out;
  }
  if (p.arithmetic != profile_.arithmetic()) {
    out.status = Status::incompatible_metadata;
    return out;
  }
  const bool needs_residual =
      p.requested & (score_bit | corrected_bit | profiled_bit | diagnostic_bit);
  const bool needs_geometry = needs_residual || (p.requested & geometry_bit);
  const bool needs_effect = needs_residual || (p.requested & effect_bit);
  const auto n = selected_.coordinates.size();
  auto bytes = p.maximum_native_bytes;
  auto fit = [&](size_t count, size_t width) {
    if (count > bytes / width)
      return false;
    bytes -= count * width;
    return true;
  };
  if (points.size() > p.maximum_models || !fit(points.size(), sizeof(Slot))) {
    out.status = Status::work_limit;
    return out;
  }
  if (n && points.size() > std::numeric_limits<size_t>::max() / n) {
    out.status = Status::work_limit;
    return out;
  }
  const auto elements = points.size() * n;
  const unsigned arrays =
      bool(p.requested & geometry_bit) + bool(p.requested & effect_bit) +
      bool(p.requested & corrected_bit) + bool(p.requested & profiled_bit);
  for (unsigned k = 0; k < arrays; ++k) {
    if (n && points.size() > bytes / sizeof(double) / n) {
      out.status = Status::work_limit;
      return out;
    }
    if (!fit(elements, sizeof(double))) {
      out.status = Status::work_limit;
      return out;
    }
  }
  if ((p.requested & diagnostic_bit) &&
      (!fit(points.size(), sizeof(statistics::ProfileResult)) ||
       !fit(elements, sizeof(size_t)))) {
    out.status = Status::work_limit;
    return out;
  }
  if (needs_geometry) {
    const auto bg = cosmology::Expansion::workspace_payload_bound(n);
    if (!bg || !fit(*bg, 1) || !fit(n, sizeof(cosmology::Request))) {
      out.status = Status::work_limit;
      return out;
    }
  }
  const bool needs_profile =
      p.requested & (score_bit | profiled_bit | diagnostic_bit);
  const size_t scratch_doubles = bool(needs_geometry) + bool(needs_effect) +
                                 bool(needs_residual) +
                                 (needs_profile ? 12 : 0);
  const size_t scratch_wide = needs_profile ? 8 : 0;
  if (!fit(n, scratch_doubles * sizeof(double) +
                  scratch_wide * sizeof(long double))) {
    out.status = Status::work_limit;
    return out;
  }
  out.status = Status::ok;
  out.slots.reserve(points.size());
  auto callbacks = p.background.maximum_callbacks,
       segments = p.background.maximum_segment_visits;
  for (const auto &point : points) {
    out.slots.emplace_back(point);
    auto &s = out.slots.back();
    auto requested_state = [](OutputState &state) {
      state.availability = cosmology::Availability::failed;
    };
    if (p.requested & geometry_bit)
      requested_state(s.geometry);
    if (p.requested & effect_bit)
      requested_state(s.effect);
    if (p.requested & corrected_bit)
      requested_state(s.corrected);
    if (p.requested & (score_bit | profiled_bit | diagnostic_bit))
      requested_state(s.profile);
    std::array<bool, 4> assessed{};
    auto assess_index = [&](OutputState *state) {
      return state == &s.geometry    ? 0u
             : state == &s.effect    ? 1u
             : state == &s.corrected ? 2u
                                     : 3u;
    };
    auto fail_pending = [&](Status status, numerics::Status cause) {
      s.status = status;
      s.numerical_status = cause;
      for (auto *state : {&s.geometry, &s.effect, &s.corrected, &s.profile})
        if (state->availability == cosmology::Availability::failed &&
            !assessed[assess_index(state)]) {
          state->status = status;
          state->numerical_status = cause;
          assessed[assess_index(state)] = true;
        }
    };
    const auto *grey = std::get_if<GreyLog1pMagnitude>(&point.source_effect);
    s.hypothesis_id = grey ? "W01/grey-log1p-zhd-magnitude/v1"
                           : "W01/no-additive-magnitude-effect/v1";
    if (grey && (!std::isfinite(grey->epsilon_mag) || grey->epsilon_mag < -.5 ||
                 grey->epsilon_mag > .5)) {
      fail_pending(Status::invalid_input,
                   std::isfinite(grey->epsilon_mag)
                       ? numerics::Status::outside_domain
                       : numerics::Status::nonfinite_input);
      continue;
    }
    auto expansion = cosmology::prepare(point.expansion, point.geometry);
    s.model_id = expansion.model_id();
    if (expansion.status() != cosmology::Status::ok) {
      s.background_status = expansion.status();
      fail_pending(Status::invalid_input,
                   expansion.status() == cosmology::Status::unsupported_domain
                       ? numerics::Status::outside_domain
                       : numerics::Status::nonfinite_input);
      continue;
    }
    auto local = p.background;
    local.maximum_callbacks = callbacks;
    local.maximum_segment_visits = segments;
    std::vector<cosmology::Request> requests;
    if (needs_geometry)
      requests.reserve(n);
    if (needs_geometry)
      for (const auto &q : selected_.coordinates)
        requests.emplace_back(
            q.z_expansion,
            (std::uint32_t)cosmology::Observable::luminosity_shape, q.observer);
    cosmology::BatchResult predicted;
    predicted.status = expansion.status();
    if (needs_geometry)
      predicted = expansion.evaluate(requests, local);
    s.background_status = predicted.status;
    s.work = predicted.work;
    if (s.work.callbacks > callbacks || s.work.segment_visits > segments) {
      fail_pending(Status::work_limit, numerics::Status::work_limit);
      continue;
    }
    callbacks -= s.work.callbacks;
    segments -= s.work.segment_visits;
    out.work.callbacks += s.work.callbacks;
    out.work.segment_visits += s.work.segment_visits;
    const bool background_valid =
        predicted.status == cosmology::Status::ok &&
        (!needs_geometry || predicted.slots.size() == n);
    if (!background_valid) {
      s.status = predicted.status == cosmology::Status::work_limit
                     ? Status::work_limit
                     : Status::numerical_failure;
      s.numerical_status = predicted.status == cosmology::Status::work_limit
                               ? numerics::Status::work_limit
                               : numerics::Status::invalid_input;
    }
    auto mark = [&](OutputState &state, bool ok, numerics::Status cause) {
      assessed[assess_index(&state)] = true;
      state.availability = ok ? cosmology::Availability::available
                              : cosmology::Availability::failed;
      state.status = ok ? Status::ok : Status::numerical_failure;
      state.numerical_status = cause;
    };
    std::vector<double> geometry, effects, base;
    bool valid = true, geometry_ok = true, effect_ok = true;
    if (needs_geometry) {
      valid = background_valid;
      geometry.reserve(n);
      if (p.requested & diagnostic_bit)
        s.background_node_indices.reserve(n);
      for (size_t i = 0; valid && i < n; ++i) {
        const auto &value = predicted.slots[i].luminosity_shape;
        if (!value.value || !(*value.value > 0)) {
          s.background_status =
              value.value ? cosmology::Status::numerical_failure : value.status;
          s.numerical_status = value.value ? numerics::Status::outside_domain
                                           : value.numerical_status;
          valid = false;
          break;
        }
        const auto shape = 5 * std::log10((long double)*value.value);
        if (!normal(shape)) {
          s.numerical_status = numerics::Status::outside_domain;
          valid = false;
          break;
        }
        geometry.push_back((double)shape);
        if (p.requested & diagnostic_bit)
          s.background_node_indices.push_back(*predicted.slots[i].node_index);
      }
      if (p.requested & geometry_bit) {
        mark(s.geometry, valid,
             valid ? numerics::Status::ok : s.numerical_status);
        if (valid)
          s.geometric_shape = geometry;
      }
      geometry_ok = valid;
      if (!valid)
        s.background_node_indices.clear();
    }
    if (needs_effect) {
      valid = true;
      effects.reserve(n);
      for (size_t i = 0; i < n; ++i) {
        const auto z = selected_.coordinates[i].z_expansion;
        const auto zmax =
            std::holds_alternative<cosmology::FixedFiveBinQ>(point.expansion)
                ? 2.5
                : 5.;
        if (z < 0 || z > zmax) {
          s.numerical_status = numerics::Status::outside_domain;
          valid = false;
          break;
        }
        double B = 0;
        if (grey && grey->epsilon_mag != 0) {
          const auto b =
              (long double)grey->epsilon_mag *
              std::log1p((long double)selected_.coordinates[i].z_expansion) /
              std::log(2.L);
          if (!normal(b)) {
            s.numerical_status = numerics::Status::outside_domain;
            valid = false;
            break;
          }
          B = (double)b;
        }
        effects.push_back(B);
      }
      if (p.requested & effect_bit) {
        mark(s.effect, valid,
             valid ? numerics::Status::ok : s.numerical_status);
        if (valid)
          s.magnitude_shifts = effects;
      }
      effect_ok = valid;
    }
    if (!geometry_ok || !effect_ok) {
      fail_pending(s.numerical_status == numerics::Status::work_limit
                       ? Status::work_limit
                       : Status::numerical_failure,
                   s.numerical_status);
      continue;
    }
    if (needs_residual) {
      base.reserve(n);
      for (size_t i = 0; i < n; ++i) {
        const auto residual = (!grey || grey->epsilon_mag == 0)
                                  ? observed_[i] - geometry[i]
                                  : observed_[i] - geometry[i] - effects[i];
        if (!std::isfinite(residual)) {
          s.numerical_status = numerics::Status::overflow;
          valid = false;
          break;
        }
        base.push_back(residual);
      }
      if (p.requested & corrected_bit)
        mark(s.corrected, valid,
             valid ? numerics::Status::ok : s.numerical_status);
      if (!valid) {
        fail_pending(Status::numerical_failure, s.numerical_status);
        continue;
      }
    }
    s.has_predictions = true;
    if (p.requested & (score_bit | profiled_bit | diagnostic_bit)) {
      auto result = profile_.evaluate(base, selected_.ordered_ids,
                                      p.maximum_forward_sensitivity);
      s.profile_status = result.status;
      mark(s.profile, result.status == statistics::DensityStatus::finite,
           result.numerical_status);
      s.numerical_status = result.numerical_status;
      if (result.status != statistics::DensityStatus::finite) {
        s.status = Status::numerical_failure;
        if (p.requested & corrected_bit)
          s.corrected_residuals = std::move(base);
        continue;
      }
      if (p.requested & score_bit)
        s.score =
            Score{result.coefficient, result.quadratic, -result.quadratic / 2};
      if (p.requested & profiled_bit)
        s.profiled_residuals = std::move(result.adjusted_residuals);
      std::vector<double>{}.swap(result.adjusted_residuals);
      if (p.requested & diagnostic_bit)
        s.diagnostics = std::move(result);
    }
    if (p.requested & corrected_bit)
      s.corrected_residuals = std::move(base);
    s.status = Status::ok;
    s.numerical_status = numerics::Status::ok;
  }
  return out;
}
} // namespace irred::supernova
