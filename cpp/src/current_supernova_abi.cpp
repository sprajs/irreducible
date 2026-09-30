#include "irred/abi.h"
#include "irred/supernova.hpp"
#include "observation_internal.hpp"
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>
namespace {
namespace sn = irred::supernova;
namespace c = irred::cosmology;
namespace num = irred::numerics;
template <class T> bool aligned(const T *p) {
  return p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0;
}
template <class T> bool descriptor(const T *p, uint64_t n, uint64_t bytes) {
  return n <= SIZE_MAX / sizeof(T) && bytes == n * sizeof(T) &&
         (!n || aligned(p));
}
bool add(size_t &bytes, size_t count, size_t width) {
  if (count > (SIZE_MAX - bytes) / width)
    return false;
  bytes += count * width;
  return true;
}
cosmo_bytes text(std::string_view s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
cosmo_strings strings(const std::vector<cosmo_bytes> &v) {
  return {v.data(), v.size(), v.size() * sizeof(cosmo_bytes)};
}
cosmo_f64_buffer doubles(const std::vector<double> &v) {
  return {sizeof(cosmo_f64_buffer), COSMO_ABI_VERSION, 2, 0, v.data(), v.size(),
          v.size() * sizeof(double)};
}
cosmo_u64_buffer indices(const std::vector<size_t> &v) {
  static_assert(sizeof(size_t) == sizeof(uint64_t));
  return {reinterpret_cast<const uint64_t *>(v.data()), v.size(),
          v.size() * sizeof(uint64_t)};
}
cosmo_output_state state(const sn::OutputState &s) {
  return {(uint32_t)s.availability, (uint32_t)s.status,
          (uint32_t)s.numerical_status, 0};
}
bool precision(uint32_t a, double b) {
  return a <= 1 && std::isfinite(b) && b > 0;
}
bool projection(const cosmo_projection_policy &p) {
  if (p.struct_size != sizeof(p) || p.abi_version != COSMO_ABI_VERSION ||
      p.has_integration > 1 || p.maximum_depth > 60 ||
      p.maximum_queries > 4096 || p.maximum_callbacks > SIZE_MAX ||
      p.maximum_segment_visits > SIZE_MAX ||
      p.maximum_evaluations_per_integral > SIZE_MAX)
    return false;
  if (!p.has_integration)
    return p.maximum_depth == 0 && p.absolute_tolerance == 0 &&
           p.relative_tolerance == 0 && p.maximum_evaluations_per_integral == 0;
  return std::isfinite(p.absolute_tolerance) &&
         std::isfinite(p.relative_tolerance) && p.absolute_tolerance >= 0 &&
         p.relative_tolerance >= 0 &&
         (p.absolute_tolerance > 0 || p.relative_tolerance > 0);
}
c::EvaluationPolicy native_projection(const cosmo_projection_policy &p,
                                      size_t bytes) {
  c::EvaluationPolicy q;
  q.maximum_queries = p.maximum_queries;
  q.maximum_callbacks = p.maximum_callbacks;
  q.maximum_segment_visits = p.maximum_segment_visits;
  q.maximum_native_bytes = bytes;
  if (p.has_integration)
    q.integration = num::IntegrationPolicy{
        p.absolute_tolerance, p.relative_tolerance,
        (size_t)p.maximum_evaluations_per_integral, p.maximum_depth};
  return q;
}
bool model_descriptor(const cosmo_current_supernova_model &m) {
  const auto &e = m.expansion;
  const auto &b = e.parameters;
  const auto &f = m.source_effect;
  const auto &a = f.parameters;
  const size_t expected = e.model == 2 ? 3 : e.model == 3 ? 5 : 1;
  return m.struct_size == sizeof(m) && m.abi_version == COSMO_ABI_VERSION &&
         !m.geometry && !m.reserved && e.struct_size == sizeof(e) &&
         e.abi_version == COSMO_ABI_VERSION && !e.reserved && e.model <= 3 &&
         b.struct_size == sizeof(b) && b.abi_version == COSMO_ABI_VERSION &&
         b.element_type == 2 && !b.reserved && b.length == expected &&
         descriptor(b.data, b.length, b.byte_length) &&
         f.struct_size == sizeof(f) && f.abi_version == COSMO_ABI_VERSION &&
         !f.reserved && f.effect <= 1 && a.struct_size == sizeof(a) &&
         a.abi_version == COSMO_ABI_VERSION && a.element_type == 2 &&
         !a.reserved && a.length == f.effect &&
         descriptor(a.data, a.length, a.byte_length);
}
c::ExpansionSpec expansion(uint32_t tag, const double *p) {
  switch (tag) {
  case 0:
    return c::LCDM(p[0]);
  case 1:
    return c::ConstantQ(p[0]);
  case 2:
    return c::CPL(p[0], p[1], p[2]);
  default:
    return c::FixedFiveBinQ({p[0], p[1], p[2], p[3], p[4]});
  }
}
struct Core {
  std::shared_ptr<const irred::observations::Prepared> source;
  sn::Consumer native;
  uint32_t arithmetic = 0;
  std::vector<cosmo_bytes> measurements, events, axes, selected_ids;
  std::vector<cosmo_magnitude_coordinate> coordinates;
  cosmo_current_supernova_view view{};
  void initialize(cosmo_observation_descriptor raw) {
    auto fill = [](const auto &input, auto &output) {
      output.reserve(input.size());
      for (const auto &s : input)
        output.push_back(text(s));
    };
    const auto &s = source->source();
    fill(s.measurement_ids, measurements);
    fill(s.event_ids, events);
    fill(s.uncertainty_axis_ids, axes);
    view.struct_size = sizeof(view);
    view.abi_version = COSMO_ABI_VERSION;
    view.status = (uint32_t)native.status();
    view.preparation_status = (uint32_t)native.preparation_status();
    view.preparation_numerical_status =
        (uint32_t)native.preparation_numerical_status();
    view.arithmetic = arithmetic;
    view.source = raw;
    view.source.measurement_ids = strings(measurements);
    view.source.event_ids = strings(events);
    view.source.uncertainty_axis_ids = strings(axes);
    const auto &sel = native.selected_source();
    fill(sel.ordered_ids, selected_ids);
    coordinates.reserve(sel.coordinates.size());
    for (const auto &q : sel.coordinates)
      coordinates.push_back(
          {sizeof(cosmo_magnitude_coordinate), COSMO_ABI_VERSION, q.z_expansion,
           q.observer.redshift, (uint32_t)q.observer.convention, 0});
    view.ordered_ids = strings(selected_ids);
    view.selected_source_indices = indices(sel.source_indices);
    view.coordinates = coordinates.data();
    view.coordinate_count = coordinates.size();
    view.coordinate_byte_length =
        coordinates.size() * sizeof(cosmo_magnitude_coordinate);
    view.arithmetic_id =
        text(arithmetic ? "F02/longdouble-cpu/v1" : "F02/binary64-legacy/v1");
    view.score_id = text(sn::Consumer::score_id);
    view.shape_convention = text(sn::Consumer::shape_convention);
    view.offset_convention = text(sn::Consumer::offset_convention);
    auto retained = native.retained_payload_bound();
    size_t bytes = sizeof(Core);
    if (!retained || !add(bytes, *retained - sizeof(native), 1) ||
        !add(bytes,
             measurements.capacity() + events.capacity() + axes.capacity() +
                 selected_ids.capacity(),
             sizeof(cosmo_bytes)) ||
        !add(bytes, coordinates.capacity(), sizeof(cosmo_magnitude_coordinate)))
      throw std::bad_alloc();
    view.retained_bytes = bytes;
  }
};
} // namespace
struct cosmo_current_supernova {
  std::shared_ptr<Core> core;
};
struct cosmo_current_supernova_result {
  std::shared_ptr<Core> core;
  sn::BatchResult native;
  std::vector<std::vector<double>> expansion_parameters, effect_parameters;
  std::vector<cosmo_current_supernova_row> rows;
  uint32_t numerical_status = (uint32_t)num::Status::ok;
};
extern "C" uint32_t cosmo_current_supernova_prepare(
    const cosmo_prepared *source, const cosmo_magnitude_selection *s,
    const cosmo_current_supernova_preparation_policy *p,
    cosmo_current_supernova **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) ||
      (!source || reinterpret_cast<uintptr_t>(source) % alignof(void *) != 0) ||
      !aligned(s) || !aligned(p))
    return COSMO_INVALID_INPUT;
  if (s->abi_version != COSMO_ABI_VERSION ||
      p->abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (s->struct_size != sizeof(*s) || p->struct_size != sizeof(*p) ||
      s->reserved || p->reserved || s->kind > 1 ||
      !precision(p->arithmetic, p->maximum_forward_sensitivity) ||
      p->maximum_selected_rows > 4096 || p->maximum_native_bytes > SIZE_MAX ||
      p->maximum_matrix_elements > SIZE_MAX ||
      p->maximum_string_bytes > SIZE_MAX ||
      !descriptor(s->source_indices.data, s->source_indices.length,
                  s->source_indices.byte_length) ||
      !descriptor(s->coordinates, s->coordinate_count,
                  s->coordinate_byte_length))
    return COSMO_INVALID_INPUT;
  if ((s->kind == 1 && (s->source_indices.length || s->coordinate_count)) ||
      (s->kind == 0 && s->source_indices.length != s->coordinate_count))
    return COSMO_INVALID_INPUT;
  try {
    auto shared = shared_native_observations(source);
    if (!shared)
      return COSMO_INVALID_INPUT;
    const auto &raw = shared->source();
    const size_t count = s->kind == 1 ? raw.values.size() : s->coordinate_count;
    if (count > 4096 ||
        (s->kind == 0 &&
         (count > p->maximum_selected_rows ||
          (count && count > p->maximum_matrix_elements / count))))
      return COSMO_INVALID_INPUT;
    auto peak = sn::preparation_payload_bound(*shared, count,
                                              (num::Arithmetic)p->arithmetic);
    size_t wrapper = sizeof(Core) + sizeof(cosmo_current_supernova);
    if (!peak || !add(wrapper, raw.values.size(), 3 * sizeof(cosmo_bytes)) ||
        !add(wrapper, count,
             sizeof(cosmo_bytes) + sizeof(cosmo_magnitude_coordinate)) ||
        wrapper > p->maximum_native_bytes ||
        *peak > p->maximum_native_bytes - wrapper)
      return COSMO_INVALID_INPUT;
    size_t chars = 0;
    for (const auto &id : raw.measurement_ids)
      if (!add(chars, id.size(), 1))
        return COSMO_INVALID_INPUT;
    for (const auto *t :
         {&raw.table_sha256, &raw.uncertainty_sha256, &raw.ordering_provenance,
          &raw.calibration_provenance, &raw.dependence_provenance})
      if (!add(chars, t->size(), 1))
        return COSMO_INVALID_INPUT;
    if (chars > p->maximum_string_bytes)
      return COSMO_INVALID_INPUT;
    sn::SelectedMagnitudeSource selected;
    if (s->kind == 1)
      selected = sn::pantheon_zhd_gt_001(shared);
    else {
      selected.source = shared;
      selected.source_indices.reserve(count);
      selected.ordered_ids.reserve(count);
      selected.coordinates.reserve(count);
      for (size_t i = 0; i < count; ++i) {
        auto index = s->source_indices.data[i];
        const auto &q = s->coordinates[i];
        if (index >= raw.values.size() ||
            (i && index <= s->source_indices.data[i - 1]) ||
            q.struct_size != sizeof(q) || q.abi_version != COSMO_ABI_VERSION ||
            q.reserved || q.observer_convention > 1)
          return COSMO_INVALID_INPUT;
        selected.source_indices.push_back(index);
        selected.ordered_ids.push_back(raw.measurement_ids[index]);
        selected.coordinates.push_back(
            {q.z_expansion,
             {q.observer_redshift, (c::Convention)q.observer_convention}});
      }
    }
    sn::PreparationPolicy policy;
    policy.arithmetic = (num::Arithmetic)p->arithmetic;
    policy.maximum_selected_rows = p->maximum_selected_rows;
    policy.maximum_matrix_elements = p->maximum_matrix_elements;
    policy.maximum_forward_sensitivity = p->maximum_forward_sensitivity;
    policy.maximum_native_bytes = p->maximum_native_bytes - wrapper;
    auto native = sn::prepare(std::move(selected), policy);
    if (native.status() == sn::Status::work_limit)
      return COSMO_INVALID_INPUT;
    cosmo_observation_descriptor raw_view{};
    if (cosmo_observation_source_view(source, &raw_view) != COSMO_OK)
      return COSMO_INVALID_INPUT;
    auto core = std::make_shared<Core>();
    core->source = std::move(shared);
    core->native = std::move(native);
    core->arithmetic = p->arithmetic;
    core->initialize(raw_view);
    auto owner = std::make_unique<cosmo_current_supernova>();
    owner->core = std::move(core);
    *out = owner.release();
    return COSMO_OK;
  } catch (const std::bad_alloc &) {
    return COSMO_ALLOCATION_FAILURE;
  } catch (...) {
    return COSMO_EXCEPTION;
  }
}
extern "C" uint32_t
cosmo_current_supernova_source_view(const cosmo_current_supernova *s,
                                    cosmo_current_supernova_view *v) {
  if (aligned(v))
    *v = {};
  if (!aligned(s) || !aligned(v))
    return COSMO_INVALID_INPUT;
  *v = s->core->view;
  return COSMO_OK;
}
extern "C" uint32_t
cosmo_current_supernova_destroy(cosmo_current_supernova *s) {
  if (!s)
    return COSMO_OK;
  if (!aligned(s))
    return COSMO_INVALID_INPUT;
  delete s;
  return COSMO_OK;
}
extern "C" uint32_t cosmo_current_supernova_evaluate(
    const cosmo_current_supernova *source,
    const cosmo_current_supernova_batch *b,
    const cosmo_current_supernova_evaluation_policy *p,
    cosmo_current_supernova_result **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(source) || !aligned(b) || !aligned(p))
    return COSMO_INVALID_INPUT;
  if (b->abi_version != COSMO_ABI_VERSION ||
      p->abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (b->struct_size != sizeof(*b) || p->struct_size != sizeof(*p) ||
      !precision(p->arithmetic, p->maximum_forward_sensitivity) ||
      !p->requested || (p->requested & ~63u) || !projection(p->projection) ||
      p->maximum_models > 65536 || p->maximum_array_elements > 1048576 ||
      p->maximum_native_bytes > SIZE_MAX || b->model_count > 65536 ||
      !descriptor(b->models, b->model_count, b->model_byte_length))
    return COSMO_INVALID_INPUT;
  const auto selected_rows = source->core->native.selected_source().coordinates.size();
  if (selected_rows && b->model_count > 1048576 / selected_rows)
    return COSMO_INVALID_INPUT;
  try {
    auto result = std::make_unique<cosmo_current_supernova_result>();
    if (source->core->native.status() != sn::Status::ok) {
      result->native.status = source->core->native.status();
      result->numerical_status =
          (uint32_t)source->core->native.preparation_numerical_status();
      result->core = source->core;
      *out = result.release();
      return COSMO_OK;
    }
    if (p->arithmetic != source->core->arithmetic) {
      result->native.status = sn::Status::incompatible_metadata;
      result->numerical_status = (uint32_t)num::Status::invalid_input;
      result->core = source->core;
      *out = result.release();
      return COSMO_OK;
    }
    const size_t count = b->model_count,
                 n = source->core->native.selected_source().coordinates.size();
    const unsigned arrays = bool(p->requested & 2u) + bool(p->requested & 4u) +
                            bool(p->requested & 8u) + bool(p->requested & 16u);
    if (n && count > SIZE_MAX / n)
      return COSMO_INVALID_INPUT;
    const auto elements = count * n;
    if (arrays && elements > SIZE_MAX / arrays)
      return COSMO_INVALID_INPUT;
    size_t bytes = sizeof(cosmo_current_supernova_result);
    if (!add(bytes, count,
             sizeof(cosmo_current_supernova_row) + sizeof(sn::ModelPoint) +
                 2 * sizeof(std::vector<double>) + 6 * sizeof(double)))
      return COSMO_INVALID_INPUT;
    if (count > p->maximum_models ||
        arrays * elements > p->maximum_array_elements ||
        bytes > p->maximum_native_bytes) {
      result->native.status = sn::Status::work_limit;
      result->numerical_status = (uint32_t)num::Status::work_limit;
      *out = result.release();
      return COSMO_OK;
    }
    for (size_t i = 0; i < count; ++i)
      if (!model_descriptor(b->models[i]))
        return COSMO_INVALID_INPUT;
    std::vector<sn::ModelPoint> models;
    models.reserve(count);
    result->expansion_parameters.reserve(count);
    result->effect_parameters.reserve(count);
    result->rows.reserve(count);
    for (size_t i = 0; i < count; ++i) {
      const auto &m = b->models[i];
      result->expansion_parameters.emplace_back(
          m.expansion.parameters.data,
          m.expansion.parameters.data + m.expansion.parameters.length);
      result->effect_parameters.emplace_back();
      if (m.source_effect.effect)
        result->effect_parameters.back().push_back(
            m.source_effect.parameters.data[0]);
      sn::SourceEffectSpec effect = sn::NoMagnitudeEffect{};
      if (m.source_effect.effect)
        effect = sn::GreyLog1pMagnitude(m.source_effect.parameters.data[0]);
      models.emplace_back(
          expansion(m.expansion.model, m.expansion.parameters.data),
          c::FlatFLRW{}, effect);
    }
    sn::Policy policy;
    policy.arithmetic = (num::Arithmetic)p->arithmetic;
    policy.maximum_forward_sensitivity = p->maximum_forward_sensitivity;
    policy.maximum_models = p->maximum_models;
    policy.maximum_native_bytes = p->maximum_native_bytes - bytes;
    policy.requested = p->requested;
    policy.background =
        native_projection(p->projection, policy.maximum_native_bytes);
    result->native = source->core->native.evaluate_batch(models, policy);
    result->core = source->core;
    if (result->native.status == sn::Status::work_limit)
      result->numerical_status = (uint32_t)num::Status::work_limit;
    else if (result->native.status != sn::Status::ok)
      result->numerical_status = (uint32_t)num::Status::invalid_input;
    for (size_t i = 0; i < result->native.slots.size(); ++i) {
      const auto &s = result->native.slots[i];
      cosmo_current_supernova_row r{};
      r.struct_size = sizeof(r);
      r.abi_version = COSMO_ABI_VERSION;
      r.model_index = i;
      r.source = b->models[i];
      r.source.expansion.parameters.data =
          result->expansion_parameters[i].data();
      r.source.source_effect.parameters.data =
          result->effect_parameters[i].data();
      r.status = (uint32_t)s.status;
      r.background_status = (uint32_t)s.background_status;
      r.numerical_status = (uint32_t)s.numerical_status;
      r.profile_status = (uint32_t)s.profile_status;
      r.geometry_state = state(s.geometry);
      r.effect_state = state(s.effect);
      r.corrected_state = state(s.corrected);
      r.profiled_state = state(s.profile);
      r.score.state = {(uint32_t)c::Availability::not_requested,
                       (uint32_t)sn::Status::invalid_input,
                       (uint32_t)num::Status::invalid_input, 0};
      r.diagnostics.state = r.score.state;
      if (p->requested & 1u) {
        r.score.state = state(s.profile);
        if (s.score) {
          r.score.offset_coefficient = s.score->offset_coefficient;
          r.score.quadratic = s.score->quadratic;
          r.score.relative_profile_score = s.score->relative_profile_score;
        }
      }
      if (p->requested & 32u) {
        r.diagnostics.state = state(s.profile);
        if (s.diagnostics) {
          const auto &d = *s.diagnostics;
          r.diagnostics.backward_residual = d.backward_residual;
          r.diagnostics.estimated_forward_sensitivity =
              d.estimated_forward_sensitivity;
          r.diagnostics.coefficient_backward_residual =
              d.coefficient_solve_backward_residual;
          r.diagnostics.coefficient_forward_sensitivity =
              d.coefficient_solve_forward_sensitivity;
          r.diagnostics.residual_l1 = d.residual_l1;
          r.diagnostics.solution_norm_inf = d.solution_norm_inf;
          r.diagnostics.adjusted_residual_l1 = d.adjusted_residual_l1;
          r.diagnostics.adjusted_solution_norm_inf =
              d.adjusted_solution_norm_inf;
        }
      }
      r.geometric_shape = doubles(s.geometric_shape);
      r.magnitude_effect = doubles(s.magnitude_shifts);
      r.corrected_residuals = doubles(s.corrected_residuals);
      r.profiled_residuals = doubles(s.profiled_residuals);
      r.background_node_indices = indices(s.background_node_indices);
      r.callbacks = s.work.callbacks;
      r.segment_visits = s.work.segment_visits;
      r.model_id = text(s.model_id);
      r.hypothesis_id = text(s.hypothesis_id);
      r.arithmetic_id = source->core->view.arithmetic_id;
      result->rows.push_back(r);
    }
    *out = result.release();
    return COSMO_OK;
  } catch (const std::bad_alloc &) {
    return COSMO_ALLOCATION_FAILURE;
  } catch (...) {
    return COSMO_EXCEPTION;
  }
}
extern "C" uint32_t cosmo_current_supernova_result_view(
    const cosmo_current_supernova_result *r,
    const cosmo_current_supernova_row **rows, uint64_t *count, uint32_t *status,
    uint32_t *numerical, uint64_t *callbacks, uint64_t *segments) {
  if (aligned(rows))
    *rows = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(status))
    *status = 0;
  if (aligned(numerical))
    *numerical = 0;
  if (aligned(callbacks))
    *callbacks = 0;
  if (aligned(segments))
    *segments = 0;
  if (!aligned(r) || !aligned(rows) || !aligned(count) || !aligned(status) ||
      !aligned(numerical) || !aligned(callbacks) || !aligned(segments))
    return COSMO_INVALID_INPUT;
  *rows = r->rows.data();
  *count = r->rows.size();
  *status = (uint32_t)r->native.status;
  *numerical = r->numerical_status;
  *callbacks = r->native.work.callbacks;
  *segments = r->native.work.segment_visits;
  return COSMO_OK;
}
extern "C" uint32_t cosmo_current_supernova_result_source_view(
    const cosmo_current_supernova_result *r, cosmo_current_supernova_view *v,
    uint32_t *available) {
  if (aligned(v))
    *v = {};
  if (aligned(available))
    *available = 0;
  if (!aligned(r) || !aligned(v) || !aligned(available))
    return COSMO_INVALID_INPUT;
  if (r->core) {
    *v = r->core->view;
    *available = 1;
  }
  return COSMO_OK;
}
extern "C" uint32_t
cosmo_current_supernova_result_destroy(cosmo_current_supernova_result *r) {
  if (!r)
    return COSMO_OK;
  if (!aligned(r))
    return COSMO_INVALID_INPUT;
  delete r;
  return COSMO_OK;
}
