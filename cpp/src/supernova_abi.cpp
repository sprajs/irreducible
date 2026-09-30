#include "irred/abi.h"
#include "irred/supernova.hpp"
#include "observation_internal.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <type_traits>
#include <vector>
namespace {
namespace sn = irred::supernova;
namespace physics = irred::cosmology;
static_assert(sizeof(cosmo_supernova_model) == 24);
static_assert(sizeof(cosmo_supernova_model_v2) == 40);
template <class T> bool aligned(const T *p) {
  return p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0;
}
static_assert(static_cast<uint32_t>(sn::Status::work_limit) ==
              COSMO_SUPERNOVA_STATUS_WORK_LIMIT);
static_assert(
    static_cast<uint32_t>(irred::numerics::Arithmetic::longdouble_cpu_v1) ==
    COSMO_SUPERNOVA_ARITHMETIC_LONGDOUBLE_CPU_V1);
cosmo_bytes bytes(std::string_view value) {
  return {reinterpret_cast<const uint8_t *>(value.data()), value.size()};
}
cosmo_f64_buffer doubles(std::span<const double> value) {
  return {sizeof(cosmo_f64_buffer),
          COSMO_ABI_VERSION,
          2,
          0,
          value.data(),
          value.size(),
          value.size_bytes()};
}
bool valid(const cosmo_supernova_policy &p) {
  return p.struct_size == sizeof(p) && p.abi_version == COSMO_ABI_VERSION &&
         !p.reserved && p.arithmetic <= 1 && p.include_residual_arrays <= 1 &&
         p.maximum_models <= 64 && p.maximum_source_rows <= 4096 &&
         p.maximum_matrix_elements <= 16777216 &&
         p.maximum_array_elements <= COSMO_MAX_BATCH_ELEMENTS &&
         p.maximum_native_output_bytes <= 512u * 1024u * 1024u &&
         p.maximum_depth > 0 && p.maximum_depth <= 64 &&
         p.maximum_evaluations_per_integral >= 3 &&
         p.maximum_total_evaluations <= std::numeric_limits<size_t>::max() &&
         p.maximum_evaluations_per_integral <=
             std::numeric_limits<size_t>::max() &&
         std::isfinite(p.absolute_tolerance) &&
         std::isfinite(p.relative_tolerance) && p.absolute_tolerance >= 0 &&
         p.relative_tolerance >= 0 &&
         (p.absolute_tolerance > 0 || p.relative_tolerance > 0) &&
         std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0;
}
sn::Policy policy(const cosmo_supernova_policy &p) {
  sn::Policy result;
  result.arithmetic = static_cast<irred::numerics::Arithmetic>(p.arithmetic);
  result.maximum_models = p.maximum_models;
  result.maximum_matrix_elements = p.maximum_matrix_elements;
  result.maximum_forward_sensitivity = p.maximum_forward_sensitivity;
  result.background.maximum_queries = p.maximum_source_rows;
  result.background.maximum_total_evaluations = p.maximum_total_evaluations;
  result.background.integration = {
      p.absolute_tolerance, p.relative_tolerance,
      static_cast<size_t>(p.maximum_evaluations_per_integral), p.maximum_depth};
  return result;
}
bool valid(const cosmo_supernova_piecewise_policy &p) {
  return p.struct_size == sizeof(p) && p.abi_version == COSMO_ABI_VERSION &&
         p.arithmetic <= 1 && p.include_residual_arrays <= 1 &&
         p.maximum_models <= 64 && p.maximum_source_rows <= 4096 &&
         p.maximum_matrix_elements <= 16777216 && p.maximum_queries <= 4096 &&
         p.maximum_array_elements <= COSMO_MAX_BATCH_ELEMENTS &&
         p.maximum_native_output_bytes <= 512u * 1024u * 1024u &&
         p.maximum_total_segment_visits <= std::numeric_limits<size_t>::max() &&
         std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0;
}
sn::Policy policy(const cosmo_supernova_piecewise_policy &p) {
  // Preparation uses only matrix/query size, arithmetic and solve quality.
  // Its unused quadrature defaults are not part of this analytic wire contract.
  sn::Policy result;
  result.arithmetic = static_cast<irred::numerics::Arithmetic>(p.arithmetic);
  result.maximum_models = p.maximum_models;
  result.maximum_matrix_elements = p.maximum_matrix_elements;
  result.maximum_forward_sensitivity = p.maximum_forward_sensitivity;
  result.background.maximum_queries = p.maximum_source_rows;
  return result;
}
bool product(size_t a, size_t b, size_t &out) {
  if (b && a > std::numeric_limits<size_t>::max() / b)
    return false;
  out = a * b;
  return true;
}
} // namespace
struct cosmo_supernova {
  sn::Consumer native;
  uint32_t arithmetic;
  std::vector<cosmo_bytes> ids;
  std::vector<uint64_t> indices;
  std::vector<double> expansion_z, observer_z;
};
struct cosmo_supernova_result_v2 {
  sn::BatchResultV2 native;
  std::vector<cosmo_supernova_slot_v2> rows;
  std::string arithmetic_id;
};
struct cosmo_supernova_result {
  sn::BatchResult native;
  std::vector<cosmo_supernova_slot> rows;
  std::string arithmetic_id;
};
template <class WirePolicy>
uint32_t prepare_supernova(const cosmo_prepared *source, const WirePolicy *p,
                           cosmo_supernova **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !source || !aligned(p))
    return COSMO_INVALID_INPUT;
  if (p->abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (!valid(*p))
    return COSMO_INVALID_INPUT;
  const auto *input = native_observations(source);
  if (!input)
    return COSMO_INVALID_INPUT;
  const auto n = input->source().values.size();
  size_t square;
  if (n > p->maximum_source_rows || !product(n, n, square) ||
      square > p->maximum_matrix_elements)
    return COSMO_INVALID_INPUT;
  try {
    auto owner = std::make_unique<cosmo_supernova>();
    owner->arithmetic = p->arithmetic;
    owner->native = sn::prepare(*input, policy(*p));
    const auto &original = owner->native.observations().source();
    const auto queries = owner->native.source_queries();
    for (auto i : owner->native.selected_source_indices()) {
      if (i >= original.measurement_ids.size() || i >= queries.size())
        return COSMO_EXCEPTION;
      owner->indices.push_back(i);
      owner->ids.push_back(bytes(original.measurement_ids[i]));
      owner->expansion_z.push_back(queries[i].z_expansion);
      owner->observer_z.push_back(queries[i].z_observer);
    }
    *out = owner.release();
    return COSMO_OK;
  } catch (const std::bad_alloc &) {
    return COSMO_ALLOCATION_FAILURE;
  } catch (...) {
    return COSMO_EXCEPTION;
  }
}
extern "C" uint32_t cosmo_supernova_prepare(const cosmo_prepared *source,
                                            const cosmo_supernova_policy *p,
                                            cosmo_supernova **out) {
  return prepare_supernova(source, p, out);
}
extern "C" uint32_t
cosmo_supernova_piecewise_prepare(const cosmo_prepared *source,
                                  const cosmo_supernova_piecewise_policy *p,
                                  cosmo_supernova **out) {
  return prepare_supernova(source, p, out);
}
extern "C" uint32_t cosmo_supernova_source_view(const cosmo_supernova *owner,
                                                cosmo_supernova_view *out) {
  if (aligned(out))
    *out = {};
  if (!aligned(owner) || !aligned(out))
    return COSMO_INVALID_INPUT;
  const auto &n = owner->native;
  const auto &m = n.probability_metadata();
  const auto &source = n.observations().source();
  out->struct_size = sizeof(*out);
  out->abi_version = COSMO_ABI_VERSION;
  out->status = static_cast<uint32_t>(n.status());
  out->preparation_status = static_cast<uint32_t>(n.preparation_status());
  out->preparation_numerical_status =
      static_cast<uint32_t>(n.preparation_numerical_status());
  out->arithmetic = owner->arithmetic;
  out->matrix_validation_assessed = n.status() == sn::Status::ok ? 1u : 0u;
  out->matrix_validation_scope =
      static_cast<uint32_t>(m.matrix_validation_scope);
  out->ordered_ids = {owner->ids.data(), owner->ids.size(),
                      owner->ids.size() * sizeof(cosmo_bytes)};
  out->selected_source_indices = {owner->indices.data(), owner->indices.size(),
                                  owner->indices.size() * sizeof(uint64_t)};
  out->source_row_count = source.values.size();
  // A failed preparation has no successfully executed factor arithmetic ID.
  if (n.status() == sn::Status::ok)
    out->arithmetic_id = bytes(m.arithmetic_id);
  out->score_id = bytes(sn::Consumer::score_id);
  out->constant_set_id = bytes(physics::Background::constants_id);
  out->radial_equation_id = bytes(physics::Background::radial_equation_id);
  out->shape_convention = bytes(sn::Consumer::shape_convention);
  out->offset_convention = bytes(sn::Consumer::offset_convention);
  out->query_provenance = bytes(n.query_provenance());
  out->table_identity = bytes(source.table_sha256);
  out->uncertainty_identity = bytes(source.uncertainty_sha256);
  out->ordering_provenance = bytes(source.ordering_provenance);
  out->calibration_provenance = bytes(source.calibration_provenance);
  out->dependence_provenance = bytes(source.dependence_provenance);
  out->source_semantics = bytes(m.source_semantics);
  out->selected_z_expansion = doubles(owner->expansion_z);
  out->selected_z_observer = doubles(owner->observer_z);
  return COSMO_OK;
}
template <class Batch, class Result, class WireModel, class WireSlot,
          class NativePoint, class NativeSlot>
uint32_t evaluate_supernova(const cosmo_supernova *owner, const Batch *b,
                            const cosmo_supernova_policy *p, Result **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(owner) || !aligned(b) || !aligned(p))
    return COSMO_INVALID_INPUT;
  if (b->abi_version != COSMO_ABI_VERSION ||
      p->abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (b->struct_size != sizeof(*b) || !valid(*p) ||
      b->model_count > p->maximum_models ||
      b->model_count > std::numeric_limits<size_t>::max() / sizeof(WireModel) ||
      b->model_byte_length != b->model_count * sizeof(WireModel) ||
      (b->model_count && (!b->models || reinterpret_cast<uintptr_t>(b->models) %
                                            alignof(WireModel))))
    return COSMO_INVALID_INPUT;
  const auto rows = owner->native.selected_source_indices().size();
  if (rows > p->maximum_source_rows)
    return COSMO_INVALID_INPUT;
  size_t elements, array_elements, array_bytes, row_bytes, temporary_bytes;
  // Native Consumer retains four ordered arrays per model (including adjusted
  // solve diagnostics), even when the borrowed export arrays are not requested.
  if (!product(b->model_count, rows, elements) ||
      !product(elements, 4, array_elements) ||
      array_elements > p->maximum_array_elements ||
      !product(array_elements, sizeof(double), array_bytes) ||
      !product(b->model_count,
               sizeof(WireSlot) + sizeof(NativeSlot) + sizeof(NativePoint),
               row_bytes) ||
      row_bytes > std::numeric_limits<size_t>::max() - array_bytes ||
      row_bytes + array_bytes > p->maximum_native_output_bytes ||
      !product(rows, sizeof(physics::Slot), temporary_bytes) ||
      temporary_bytes >
          std::numeric_limits<size_t>::max() - row_bytes - array_bytes)
    return COSMO_INVALID_INPUT;
  for (uint64_t i = 0; i < b->model_count; ++i)
    if (b->models[i].reserved)
      return COSMO_INVALID_INPUT;
  try {
    std::vector<NativePoint> points;
    points.reserve(b->model_count);
    for (uint64_t i = 0; i < b->model_count; ++i) {
      const auto &m = b->models[i];
      if constexpr (std::is_same_v<NativePoint, sn::ModelPointV2>)
        points.emplace_back(static_cast<physics::Model>(m.model), m.omega_m,
                            m.constant_q, m.w0, m.wa);
      else
        points.push_back(
            {static_cast<physics::Model>(m.model <= 1 ? m.model : UINT32_MAX),
             m.omega_m, m.constant_q});
    }
    auto result = std::make_unique<Result>();
    if constexpr (std::is_same_v<NativePoint, sn::ModelPointV2>)
      result->native = owner->native.evaluate_batch_v2(points, policy(*p));
    else
      result->native = owner->native.evaluate_batch(points, policy(*p));
    if (owner->native.status() == sn::Status::ok)
      result->arithmetic_id =
          owner->native.probability_metadata().arithmetic_id;
    if (result->native.slots.size() > b->model_count)
      return COSMO_EXCEPTION;
    result->rows.reserve(result->native.slots.size());
    uint64_t used = 0;
    for (size_t i = 0; i < result->native.slots.size(); ++i) {
      const auto &s = [&]() -> const sn::Slot & {
        if constexpr (std::is_same_v<NativePoint, sn::ModelPointV2>)
          return result->native.slots[i].calculation;
        else
          return result->native.slots[i];
      }();
      const auto &d = s.solve_diagnostics;
      if (s.background_evaluations > p->maximum_total_evaluations - used)
        return COSMO_EXCEPTION;
      used += s.background_evaluations;
      WireSlot row{};
      row.model_index = i;
      row.source_parameters = b->models[i];
      row.status = static_cast<uint32_t>(s.status);
      row.background_status = static_cast<uint32_t>(s.background_status);
      row.numerical_status = static_cast<uint32_t>(s.numerical_status);
      row.profile_status = static_cast<uint32_t>(s.profile_status);
      row.background_evaluations = s.background_evaluations;
      auto model = [&] {
        if constexpr (std::is_same_v<NativePoint, sn::ModelPointV2>) {
          const auto &m = points[i];
          if (m.model == physics::Model::flat_cpl_late_v1 && m.constant_q == 0)
            return physics::prepare_cpl(
                {sn::Consumer::computational_h0_km_s_mpc, m.omega_m, m.w0,
                 m.wa});
          auto value = physics::Parameters{
              m.model, sn::Consumer::computational_h0_km_s_mpc, m.omega_m,
              m.constant_q};
          value.w0 = m.w0;
          value.wa = m.wa;
          return physics::prepare(value);
        } else
          return physics::prepare({s.source.model,
                                   sn::Consumer::computational_h0_km_s_mpc,
                                   s.source.omega_m, s.source.constant_q});
      }();
      row.model_id = bytes(model.model_id());
      row.radial_equation_id = bytes(physics::Background::radial_equation_id);
      row.score_id = bytes(sn::Consumer::score_id);
      if (owner->native.status() == sn::Status::ok)
        row.arithmetic_id = bytes(result->arithmetic_id);
      if (s.status == sn::Status::ok) {
        row.offset_coefficient = s.offset_coefficient;
        row.quadratic = s.quadratic;
        row.relative_profile_score = s.relative_profile_score;
        row.backward_residual = d.backward_residual;
        row.estimated_forward_sensitivity = d.estimated_forward_sensitivity;
        row.coefficient_solve_backward_residual =
            d.coefficient_solve_backward_residual;
        row.coefficient_solve_forward_sensitivity =
            d.coefficient_solve_forward_sensitivity;
        row.residual_l1 = d.residual_l1;
        row.solution_norm_inf = d.solution_norm_inf;
        row.adjusted_residual_l1 = d.adjusted_residual_l1;
        row.adjusted_solution_norm_inf = d.adjusted_solution_norm_inf;
        if (p->include_residual_arrays) {
          row.shape_magnitudes = doubles(s.shape_magnitudes);
          row.base_residuals = doubles(s.base_residuals);
          row.profiled_residuals = doubles(s.profiled_residuals);
        }
      }
      result->rows.push_back(row);
    }
    *out = result.release();
    return COSMO_OK;
  } catch (const std::bad_alloc &) {
    return COSMO_ALLOCATION_FAILURE;
  } catch (...) {
    return COSMO_EXCEPTION;
  }
}
extern "C" uint32_t cosmo_supernova_evaluate(const cosmo_supernova *owner,
                                             const cosmo_supernova_batch *b,
                                             const cosmo_supernova_policy *p,
                                             cosmo_supernova_result **out) {
  return evaluate_supernova<cosmo_supernova_batch, cosmo_supernova_result,
                            cosmo_supernova_model, cosmo_supernova_slot,
                            sn::ModelPoint, sn::Slot>(owner, b, p, out);
}
extern "C" uint32_t cosmo_supernova_evaluate_v2(
    const cosmo_supernova *owner, const cosmo_supernova_batch_v2 *b,
    const cosmo_supernova_policy *p, cosmo_supernova_result_v2 **out) {
  return evaluate_supernova<cosmo_supernova_batch_v2, cosmo_supernova_result_v2,
                            cosmo_supernova_model_v2, cosmo_supernova_slot_v2,
                            sn::ModelPointV2, sn::SlotV2>(owner, b, p, out);
}
extern "C" uint32_t
cosmo_supernova_result_view(const cosmo_supernova_result *owner,
                            const cosmo_supernova_slot **rows, uint64_t *count,
                            uint32_t *status) {
  if (aligned(rows))
    *rows = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(status))
    *status = COSMO_SUPERNOVA_STATUS_INVALID_INPUT;
  if (!aligned(owner) || !aligned(rows) || !aligned(count) || !aligned(status))
    return COSMO_INVALID_INPUT;
  *rows = owner->rows.data();
  *count = owner->rows.size();
  *status = static_cast<uint32_t>(owner->native.status);
  return COSMO_OK;
}
extern "C" uint32_t cosmo_supernova_destroy(cosmo_supernova *owner) {
  delete owner;
  return COSMO_OK;
}
extern "C" uint32_t
cosmo_supernova_result_destroy(cosmo_supernova_result *owner) {
  delete owner;
  return COSMO_OK;
}

extern "C" uint32_t
cosmo_supernova_result_v2_view(const cosmo_supernova_result_v2 *owner,
                               const cosmo_supernova_slot_v2 **rows,
                               uint64_t *count, uint32_t *status) {
  if (aligned(rows))
    *rows = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(status))
    *status = COSMO_SUPERNOVA_STATUS_INVALID_INPUT;
  if (!aligned(owner) || !aligned(rows) || !aligned(count) || !aligned(status))
    return COSMO_INVALID_INPUT;
  *rows = owner->rows.data();
  *count = owner->rows.size();
  *status = static_cast<uint32_t>(owner->native.status);
  return COSMO_OK;
}
extern "C" uint32_t
cosmo_supernova_result_v2_destroy(cosmo_supernova_result_v2 *owner) {
  delete owner;
  return COSMO_OK;
}

struct cosmo_supernova_piecewise_result {
  sn::PiecewiseBatchResult native;
  std::vector<cosmo_supernova_piecewise_slot> rows;
  std::string arithmetic_id;
  std::vector<std::string> id_storage;
  std::vector<cosmo_bytes> ids;
  std::vector<uint64_t> indices;
  std::vector<double> expansion_z, observer_z;
};

extern "C" uint32_t
cosmo_supernova_piecewise_evaluate(const cosmo_supernova *owner,
                                   const cosmo_supernova_piecewise_batch *b,
                                   const cosmo_supernova_piecewise_policy *p,
                                   cosmo_supernova_piecewise_result **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(owner) || !aligned(b) || !aligned(p))
    return COSMO_INVALID_INPUT;
  if (b->abi_version != COSMO_ABI_VERSION ||
      p->abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (b->struct_size != sizeof(*b) || !valid(*p) || b->model_count > 64 ||
      b->model_count > std::numeric_limits<size_t>::max() /
                           sizeof(cosmo_supernova_piecewise_model) ||
      b->model_bytes !=
          b->model_count * sizeof(cosmo_supernova_piecewise_model) ||
      (b->model_count && !aligned(b->models)))
    return COSMO_INVALID_INPUT;
  auto exhausted = [&]() -> uint32_t {
    try {
      auto result = std::make_unique<cosmo_supernova_piecewise_result>();
      result->native.status = sn::Status::work_limit;
      *out = result.release();
      return COSMO_OK;
    } catch (const std::bad_alloc &) {
      return COSMO_ALLOCATION_FAILURE;
    } catch (...) {
      return COSMO_EXCEPTION;
    }
  };
  if (b->model_count > p->maximum_models)
    return exhausted();
  const auto n = owner->indices.size();
  size_t total = sizeof(cosmo_supernova_piecewise_result);
  auto add = [&](size_t count, size_t width) {
    size_t value;
    if (!product(count, width, value) ||
        value > std::numeric_limits<size_t>::max() - total)
      return false;
    total += value;
    return true;
  };
  size_t elements, arrays;
  // Retained native arrays exist even when their borrowed export is disabled.
  // Include metadata copies, point/row storage and peak linear native
  // workspace. Caller buffers and allocator bookkeeping are outside this
  // aggregate cap.
  if (!product(b->model_count, n, elements) || !product(elements, 4, arrays) ||
      !add(arrays, sizeof(double)) ||
      !add(b->model_count, sizeof(cosmo_supernova_piecewise_slot) +
                               sizeof(sn::PiecewiseSlot) +
                               sizeof(sn::PiecewiseModelPoint)) ||
      !add(n, sizeof(std::string) + sizeof(cosmo_bytes) + sizeof(uint64_t) +
                  2 * sizeof(double)) ||
      !add(n, sizeof(physics::PiecewiseSlot) + sizeof(physics::Query) +
                  12 * sizeof(double) + 4 * sizeof(long double)))
    return COSMO_INVALID_INPUT;
  for (const auto &id : owner->ids)
    if (!add(id.length, 1) || !add(1, 1))
      return COSMO_INVALID_INPUT;
  if (!add(owner->native.probability_metadata().arithmetic_id.size(), 1) ||
      !add(1, 1))
    return COSMO_INVALID_INPUT;
  if (arrays > p->maximum_array_elements ||
      total > p->maximum_native_output_bytes)
    return exhausted();
  try {
    std::vector<sn::PiecewiseModelPoint> points;
    points.reserve(b->model_count);
    for (size_t i = 0; i < b->model_count; ++i) {
      std::array<double, 5> q;
      std::copy_n(b->models[i].q, 5, q.begin());
      points.emplace_back(q);
    }
    auto result = std::make_unique<cosmo_supernova_piecewise_result>();
    result->indices = owner->indices;
    result->expansion_z = owner->expansion_z;
    result->observer_z = owner->observer_z;
    result->id_storage.reserve(n);
    for (const auto &id : owner->ids)
      result->id_storage.emplace_back(reinterpret_cast<const char *>(id.data),
                                      id.length);
    result->ids.reserve(n);
    for (const auto &id : result->id_storage)
      result->ids.push_back(bytes(id));
    if (owner->native.status() == sn::Status::ok)
      result->arithmetic_id =
          owner->native.probability_metadata().arithmetic_id;
    sn::PiecewiseEvaluationPolicy native_policy;
    native_policy.arithmetic =
        static_cast<irred::numerics::Arithmetic>(p->arithmetic);
    native_policy.maximum_models = p->maximum_models;
    native_policy.maximum_forward_sensitivity = p->maximum_forward_sensitivity;
    native_policy.background.maximum_queries = p->maximum_queries;
    native_policy.background.maximum_segment_visits =
        p->maximum_total_segment_visits;
    result->native =
        owner->native.evaluate_piecewise_batch(points, native_policy);
    if (result->native.slots.size() > b->model_count)
      return COSMO_EXCEPTION;
    result->rows.reserve(result->native.slots.size());
    size_t used = 0;
    for (size_t i = 0; i < result->native.slots.size(); ++i) {
      const auto &s = result->native.slots[i];
      const auto &d = s.solve_diagnostics;
      if (s.segment_visits > p->maximum_total_segment_visits - used)
        return COSMO_EXCEPTION;
      used += s.segment_visits;
      cosmo_supernova_piecewise_slot row{};
      row.struct_size = sizeof(row);
      row.abi_version = COSMO_ABI_VERSION;
      row.model_index = i;
      row.source_parameters = b->models[i];
      row.status = static_cast<uint32_t>(s.status);
      row.background_status = static_cast<uint32_t>(s.background_status);
      row.numerical_status = static_cast<uint32_t>(s.numerical_status);
      row.profile_status = static_cast<uint32_t>(s.profile_status);
      row.segment_visits = s.segment_visits;
      row.model_id = bytes(physics::PiecewiseBackground::model_id);
      row.radial_equation_id = bytes(physics::Background::radial_equation_id);
      row.score_id = bytes(sn::Consumer::score_id);
      row.arithmetic_id = bytes(result->arithmetic_id);
      row.ordered_ids = {result->ids.data(), n, n * sizeof(cosmo_bytes)};
      row.selected_source_indices = result->indices.data();
      row.selected_count = n;
      row.selected_index_bytes = n * sizeof(uint64_t);
      row.expansion_z = doubles(result->expansion_z);
      row.observer_z = doubles(result->observer_z);
      if (s.status == sn::Status::ok) {
        row.has_profile_payload = 1;
        row.offset_coefficient = s.offset_coefficient;
        row.quadratic = s.quadratic;
        row.relative_profile_score = s.relative_profile_score;
        row.backward_residual = d.backward_residual;
        row.estimated_forward_sensitivity = d.estimated_forward_sensitivity;
        row.coefficient_solve_backward_residual =
            d.coefficient_solve_backward_residual;
        row.coefficient_solve_forward_sensitivity =
            d.coefficient_solve_forward_sensitivity;
        row.residual_l1 = d.residual_l1;
        row.solution_norm_inf = d.solution_norm_inf;
        row.adjusted_residual_l1 = d.adjusted_residual_l1;
        row.adjusted_solution_norm_inf = d.adjusted_solution_norm_inf;
        if (p->include_residual_arrays) {
          row.shape_magnitudes = doubles(s.shape_magnitudes);
          row.base_residuals = doubles(s.base_residuals);
          row.profiled_residuals = doubles(s.profiled_residuals);
        }
      }
      result->rows.push_back(row);
    }
    if (used != result->native.segment_visits)
      return COSMO_EXCEPTION;
    *out = result.release();
    return COSMO_OK;
  } catch (const std::bad_alloc &) {
    return COSMO_ALLOCATION_FAILURE;
  } catch (...) {
    return COSMO_EXCEPTION;
  }
}
extern "C" uint32_t cosmo_supernova_piecewise_result_view(
    const cosmo_supernova_piecewise_result *owner,
    const cosmo_supernova_piecewise_slot **rows, uint64_t *count,
    uint32_t *status, uint64_t *segments) {
  if (aligned(rows))
    *rows = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(status))
    *status = COSMO_SUPERNOVA_STATUS_INVALID_INPUT;
  if (aligned(segments))
    *segments = 0;
  if (!aligned(owner) || !aligned(rows) || !aligned(count) ||
      !aligned(status) || !aligned(segments))
    return COSMO_INVALID_INPUT;
  *rows = owner->rows.data();
  *count = owner->rows.size();
  *status = static_cast<uint32_t>(owner->native.status);
  *segments = owner->native.segment_visits;
  return COSMO_OK;
}
extern "C" uint32_t cosmo_supernova_piecewise_result_destroy(
    cosmo_supernova_piecewise_result *owner) {
  delete owner;
  return COSMO_OK;
}
