#include "irred/abi.h"
#include "irred/bao.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <vector>
namespace {
namespace b = irred::bao;
namespace c = irred::cosmology;
namespace n = irred::numerics;
static_assert(static_cast<uint32_t>(b::Observable::transverse_over_ruler) ==
              COSMO_BAO_OBSERVABLE_DM_OVER_RS);
static_assert(static_cast<uint32_t>(b::Observable::hubble_over_ruler) ==
              COSMO_BAO_OBSERVABLE_DH_OVER_RS);
static_assert(static_cast<uint32_t>(b::Observable::volume_over_ruler) ==
              COSMO_BAO_OBSERVABLE_DV_OVER_RS);
static_assert(static_cast<uint32_t>(b::RowRole::synthetic_control) ==
              COSMO_BAO_ROLE_SYNTHETIC_CONTROL);
static_assert(static_cast<uint32_t>(c::Model::flat_cpl_late_v1) ==
              COSMO_BACKGROUND_V2_MODEL_FLAT_CPL_LATE_V1);
constexpr std::string_view redshift_id = "P01/released-effective-redshift/v1";
constexpr std::string_view h0_id = "P01/computational-H0-fixed-70-km-s-Mpc/v1";
template <class T> bool aligned(const T *p) {
  return p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0;
}
bool multiply(uint64_t a, uint64_t v, uint64_t &result) {
  if (v && a > UINT64_MAX / v)
    return false;
  result = a * v;
  return true;
}
bool consume(uint64_t count, uint64_t size, uint64_t &left) {
  uint64_t bytes;
  if (!multiply(count, size, bytes) || bytes > left)
    return false;
  left -= bytes;
  return true;
}
template <class T>
bool bounded(const T *p, uint64_t count, uint64_t bytes, uint64_t cap) {
  uint64_t expected;
  return count <= cap && count <= SIZE_MAX / sizeof(T) &&
         multiply(count, sizeof(T), expected) && bytes == expected &&
         (!count || aligned(p));
}
bool valid(const cosmo_bao_policy &p) {
  const auto &q = p.background;
  return p.struct_size == sizeof(p) && p.abi_version == COSMO_ABI_VERSION &&
         !p.reserved && p.arithmetic <= 1 && p.include_predictions <= 1 &&
         p.include_residuals <= 1 && p.maximum_models <= 64 &&
         p.maximum_rows <= 4096 && p.maximum_matrix_elements <= 16777216 &&
         p.maximum_string_bytes <= 16777216 &&
         p.maximum_native_bytes <= 512u * 1024u * 1024u &&
         p.maximum_array_elements <= COSMO_MAX_BATCH_ELEMENTS &&
         p.maximum_native_output_bytes <= 512u * 1024u * 1024u &&
         std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0 && q.struct_size == sizeof(q) &&
         q.abi_version == COSMO_ABI_VERSION && !q.reserved &&
         q.maximum_depth > 0 && q.maximum_depth <= 64 &&
         q.maximum_evaluations_per_integral >= 3 &&
         q.maximum_evaluations_per_integral <= SIZE_MAX &&
         q.maximum_total_evaluations <= SIZE_MAX &&
         q.maximum_parameters <= 64 && q.maximum_queries <= 4096 &&
         q.maximum_slots <= COSMO_MAX_BATCH_ELEMENTS &&
         q.maximum_native_output_bytes <= 512u * 1024u * 1024u &&
         std::isfinite(q.absolute_tolerance) &&
         std::isfinite(q.relative_tolerance) && q.absolute_tolerance >= 0 &&
         q.relative_tolerance >= 0 &&
         (q.absolute_tolerance > 0 || q.relative_tolerance > 0);
}
b::DensityPolicy policy(const cosmo_bao_policy &p) {
  b::DensityPolicy q{};
  q.maximum_models = p.maximum_models;
  q.maximum_matrix_elements = p.maximum_matrix_elements;
  q.maximum_string_bytes = p.maximum_string_bytes;
  q.maximum_native_bytes = p.maximum_native_bytes;
  q.maximum_forward_sensitivity = p.maximum_forward_sensitivity;
  q.arithmetic = static_cast<n::Arithmetic>(p.arithmetic);
  q.observables.maximum_queries =
      std::min(p.maximum_rows, p.background.maximum_queries);
  q.observables.maximum_native_bytes = p.background.maximum_native_output_bytes;
  q.observables.background.maximum_queries = q.observables.maximum_queries;
  q.observables.background.maximum_total_evaluations =
      p.background.maximum_total_evaluations;
  q.observables.background.integration = {
      p.background.absolute_tolerance, p.background.relative_tolerance,
      static_cast<size_t>(p.background.maximum_evaluations_per_integral),
      p.background.maximum_depth};
  return q;
}
cosmo_bytes view(std::string_view s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
cosmo_f64_buffer view(std::span<const double> s) {
  return {sizeof(cosmo_f64_buffer),
          COSMO_ABI_VERSION,
          2,
          0,
          s.data(),
          s.size(),
          s.size_bytes()};
}
bool doubles(const cosmo_f64_buffer &s, uint64_t count) {
  return s.struct_size == sizeof(s) && s.abi_version == COSMO_ABI_VERSION &&
         s.element_type == 2 && !s.reserved && s.length == count &&
         bounded(s.data, s.length, s.byte_length, count);
}
bool string(const cosmo_bytes &s, uint64_t &left) {
  return bounded(s.data, s.length, s.length, left) &&
         consume(s.length, 1, left);
}
bool strings(const cosmo_strings &s, uint64_t count, uint64_t &left) {
  if (s.length != count || !bounded(s.data, s.length, s.byte_length, count))
    return false;
  for (uint64_t i = 0; i < count; ++i)
    if (!string(s.data[i], left))
      return false;
  return true;
}
std::string copy(cosmo_bytes s) {
  return s.length
             ? std::string(reinterpret_cast<const char *>(s.data), s.length)
             : std::string{};
}
std::string_view text(cosmo_bytes s) {
  return s.length ? std::string_view(reinterpret_cast<const char *>(s.data),
                                     s.length)
                  : std::string_view{};
}
} // namespace
struct cosmo_bao {
  b::DensityInput source;
  b::PreparedDensity native;
  cosmo_bao_policy preparation{};
  std::vector<cosmo_bao_query> queries;
  std::vector<cosmo_bytes> ids;
  std::string redshift, ruler, h0;
};
struct cosmo_bao_result {
  b::DensityBatch native;
  std::vector<cosmo_bao_row> rows;
  std::vector<std::string> ordered_ids;
  std::vector<cosmo_bytes> id_views;
  std::vector<cosmo_bao_query> queries;
};
extern "C" uint32_t cosmo_bao_prepare(const cosmo_bao_descriptor *d,
                                      const cosmo_bao_policy *p,
                                      cosmo_bao **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(d) || !aligned(p))
    return COSMO_INVALID_INPUT;
  if (d->abi_version != COSMO_ABI_VERSION ||
      p->abi_version != COSMO_ABI_VERSION ||
      p->background.abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (d->struct_size != sizeof(*d) || !valid(*p) || d->role > 1 ||
      d->covariance_unit != 0)
    return COSMO_INVALID_INPUT;
  const auto count = d->query_count;
  uint64_t square;
  if (!count || count > p->maximum_rows || !multiply(count, count, square) ||
      square > p->maximum_matrix_elements ||
      !bounded(d->queries, count, d->query_byte_length, p->maximum_rows) ||
      !doubles(d->observed, count) || !doubles(d->covariance, square))
    return COSMO_INVALID_INPUT;
  uint64_t strings_left = p->maximum_string_bytes;
  for (const auto *s : {&d->table_identity, &d->covariance_identity,
                        &d->ordering_provenance, &d->calibration_provenance,
                        &d->dependence_provenance, &d->redshift_convention,
                        &d->ruler_convention, &d->computational_h0_convention})
    if (!string(*s, strings_left))
      return COSMO_INVALID_INPUT;
  if (!strings(d->ordered_ids, count, strings_left) ||
      !strings(d->covariance_axis_ids, count, strings_left) ||
      text(d->redshift_convention) != redshift_id ||
      text(d->ruler_convention) != b::ruler_convention_id ||
      text(d->computational_h0_convention) != h0_id)
    return COSMO_INVALID_INPUT;
  for (uint64_t i = 0; i < count; ++i)
    if (d->queries[i].reserved ||
        text(d->ordered_ids.data[i]) != text(d->covariance_axis_ids.data[i]))
      return COSMO_INVALID_INPUT;
  // Bound wrapper retained inputs plus native copies/factors/workspace before
  // copying.
  uint64_t left = p->maximum_native_bytes;
  if (!consume(1, sizeof(cosmo_bao), left) ||
      !consume(count,
               sizeof(cosmo_bao_query) + 2 * sizeof(b::Query) +
                   4 * sizeof(std::string) + 2 * sizeof(cosmo_bytes) +
                   12 * sizeof(long double) + 8 * sizeof(double),
               left) ||
      !consume(square, 6 * sizeof(double) + sizeof(long double), left) ||
      !consume(p->maximum_string_bytes - strings_left, 3, left))
    return COSMO_INVALID_INPUT;
  try {
    auto owner = std::make_unique<cosmo_bao>();
    owner->preparation = *p;
    owner->queries.assign(d->queries, d->queries + count);
    auto &s = owner->source;
    s.role = static_cast<b::RowRole>(d->role);
    s.covariance_unit = static_cast<b::CovarianceUnit>(d->covariance_unit);
    s.table_identity = copy(d->table_identity);
    s.covariance_identity = copy(d->covariance_identity);
    s.ordering_provenance = copy(d->ordering_provenance);
    s.calibration_provenance = copy(d->calibration_provenance);
    s.dependence_provenance = copy(d->dependence_provenance);
    owner->redshift = copy(d->redshift_convention);
    owner->ruler = copy(d->ruler_convention);
    owner->h0 = copy(d->computational_h0_convention);
    s.queries.reserve(count);
    s.ordered_ids.reserve(count);
    owner->ids.reserve(count);
    for (uint64_t i = 0; i < count; ++i) {
      s.queries.push_back({d->queries[i].z, static_cast<b::Observable>(
                                                d->queries[i].observable)});
      s.ordered_ids.push_back(copy(d->ordered_ids.data[i]));
    }
    s.observed.assign(d->observed.data, d->observed.data + count);
    s.covariance.assign(d->covariance.data, d->covariance.data + square);
    owner->native = b::prepare_density(s, policy(*p));
    for (const auto &id : s.ordered_ids)
      owner->ids.push_back(view(id));
    *out = owner.release();
    return COSMO_OK;
  } catch (const std::bad_alloc &) {
    return COSMO_ALLOCATION_FAILURE;
  } catch (...) {
    return COSMO_EXCEPTION;
  }
}
extern "C" uint32_t cosmo_bao_source_view(const cosmo_bao *owner,
                                          cosmo_bao_source_view_t *out) {
  if (aligned(out))
    *out = {};
  if (!owner || !aligned(out))
    return COSMO_INVALID_INPUT;
  const auto &s = owner->source;
  out->struct_size = sizeof(*out);
  out->abi_version = COSMO_ABI_VERSION;
  out->status = static_cast<uint32_t>(owner->native.status());
  out->numerical_status =
      static_cast<uint32_t>(owner->native.numerical_status());
  out->prepare_policy = owner->preparation;
  out->arithmetic_id = view(owner->native.metadata().arithmetic_id);
  out->equation_id = view(b::equation_id);
  out->constants_id = view(irred::constant_set_id);
  auto &d = out->source;
  d.struct_size = sizeof(d);
  d.abi_version = COSMO_ABI_VERSION;
  d.role = static_cast<uint32_t>(s.role);
  d.covariance_unit = static_cast<uint32_t>(s.covariance_unit);
  d.queries = owner->queries.data();
  d.query_count = owner->queries.size();
  d.query_byte_length = owner->queries.size() * sizeof(cosmo_bao_query);
  d.observed = view(s.observed);
  d.covariance = view(s.covariance);
  d.ordered_ids = {owner->ids.data(), owner->ids.size(),
                   owner->ids.size() * sizeof(cosmo_bytes)};
  d.covariance_axis_ids = d.ordered_ids;
  d.table_identity = view(s.table_identity);
  d.covariance_identity = view(s.covariance_identity);
  d.ordering_provenance = view(s.ordering_provenance);
  d.calibration_provenance = view(s.calibration_provenance);
  d.dependence_provenance = view(s.dependence_provenance);
  d.redshift_convention = view(owner->redshift);
  d.ruler_convention = view(owner->ruler);
  d.computational_h0_convention = view(owner->h0);
  return COSMO_OK;
}
extern "C" uint32_t cosmo_bao_destroy(cosmo_bao *owner) {
  delete owner;
  return COSMO_OK;
}
extern "C" uint32_t cosmo_bao_evaluate(const cosmo_bao *owner,
                                       const cosmo_bao_batch *batch,
                                       const cosmo_bao_policy *p,
                                       cosmo_bao_result **out) {
  if (aligned(out))
    *out = nullptr;
  if (!owner || !aligned(batch) || !aligned(p) || !aligned(out))
    return COSMO_INVALID_INPUT;
  if (batch->abi_version != COSMO_ABI_VERSION ||
      p->abi_version != COSMO_ABI_VERSION ||
      p->background.abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (batch->struct_size != sizeof(*batch) || !valid(*p) ||
      !bounded(batch->models, batch->model_count, batch->model_byte_length,
               p->maximum_models) ||
      batch->model_count > p->background.maximum_parameters)
    return COSMO_INVALID_INPUT;
  const auto count = batch->model_count, rows = owner->source.queries.size();
  uint64_t arrays;
  uint64_t left = p->maximum_native_output_bytes;
  if (!multiply(count, rows, arrays) ||
      arrays > p->maximum_array_elements / 2 ||
      arrays > p->background.maximum_slots ||
      !consume(1, sizeof(cosmo_bao_result), left) ||
      !consume(count,
               sizeof(cosmo_bao_row) + sizeof(b::DensitySlot) +
                   sizeof(b::ModelQuery),
               left) ||
      !consume(count, 256, left) ||
      !consume(rows,
               sizeof(std::string) + sizeof(cosmo_bytes) +
                   sizeof(cosmo_bao_query),
               left) ||
      !consume(arrays, 2 * sizeof(double), left))
    return COSMO_INVALID_INPUT;
  for (const auto &id : owner->source.ordered_ids)
    if (!consume(id.size() + 1, 1, left))
      return COSMO_INVALID_INPUT;
  for (uint64_t i = 0; i < count; ++i)
    if (batch->models[i].parameters.reserved)
      return COSMO_INVALID_INPUT;
  try {
    auto result = std::make_unique<cosmo_bao_result>();
    result->ordered_ids = owner->source.ordered_ids;
    result->queries = owner->queries;
    result->id_views.reserve(result->ordered_ids.size());
    for (const auto &id : result->ordered_ids)
      result->id_views.push_back(view(id));
    std::vector<b::ModelQuery> models;
    models.reserve(count);
    for (uint64_t i = 0; i < count; ++i) {
      const auto &v = batch->models[i];
      const auto &m = v.parameters;
      const auto model = static_cast<c::Model>(m.model);
      auto prepared =
          (m.model == COSMO_BACKGROUND_V2_MODEL_FLAT_CPL_LATE_V1 &&
           m.constant_q == 0)
              ? c::prepare_cpl(c::CplParameters{70, m.omega_m, m.w0, m.wa})
              : c::prepare(c::Parameters{model, 70, m.omega_m, m.constant_q,
                                         m.w0, m.wa});
      models.push_back({std::move(prepared), b::Ruler(v.h0_rd_km_s)});
    }
    result->native = owner->native.evaluate(models, policy(*p));
    result->rows.resize(result->native.slots.size());
    for (size_t i = 0; i < result->rows.size(); ++i) {
      auto &row = result->rows[i];
      const auto &slot = result->native.slots[i];
      row.model_index = i;
      row.ordered_ids = {result->id_views.data(), result->id_views.size(),
                         result->id_views.size() * sizeof(cosmo_bytes)};
      row.queries = result->queries.data();
      row.query_count = result->queries.size();
      row.query_byte_length = result->queries.size() * sizeof(cosmo_bao_query);
      row.source_parameters = batch->models[i];
      row.status = static_cast<uint32_t>(slot.result.density.status);
      row.background_status = static_cast<uint32_t>(slot.background_status);
      row.numerical_status = static_cast<uint32_t>(slot.numerical_status);
      row.model_id = view(slot.model_id);
      row.arithmetic_id = view(slot.arithmetic_id);
      row.equation_id = view(slot.equation_id);
      row.ruler_convention_id = view(b::ruler_convention_id);
      row.computational_h0_convention_id = view(h0_id);
      row.evaluations = slot.evaluations;
      if (slot.result.density.status ==
          irred::statistics::DensityStatus::finite) {
        row.log_density = slot.result.density.log_value;
        row.quadratic = slot.result.quadratic;
        row.log_determinant = slot.result.log_determinant;
        row.normalization = slot.result.normalization;
        row.backward_residual = slot.result.backward_residual;
        row.estimated_forward_sensitivity =
            slot.result.estimated_forward_sensitivity;
        if (p->include_predictions)
          row.predictions = view(slot.predictions);
        if (p->include_residuals)
          row.residuals = view(slot.residuals);
      }
    }
    *out = result.release();
    return COSMO_OK;
  } catch (const std::bad_alloc &) {
    return COSMO_ALLOCATION_FAILURE;
  } catch (...) {
    return COSMO_EXCEPTION;
  }
}
extern "C" uint32_t cosmo_bao_result_view(const cosmo_bao_result *owner,
                                          const cosmo_bao_row **rows,
                                          uint64_t *count, uint32_t *status,
                                          uint32_t *numerical,
                                          uint64_t *evaluations) {
  if (aligned(rows))
    *rows = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(status))
    *status =
        static_cast<uint32_t>(irred::statistics::DensityStatus::invalid_input);
  if (aligned(numerical))
    *numerical = static_cast<uint32_t>(n::Status::invalid_input);
  if (aligned(evaluations))
    *evaluations = 0;
  if (!owner || !aligned(rows) || !aligned(count) || !aligned(status) ||
      !aligned(numerical) || !aligned(evaluations))
    return COSMO_INVALID_INPUT;
  *rows = owner->rows.data();
  *count = owner->rows.size();
  *status = static_cast<uint32_t>(owner->native.status);
  *numerical = static_cast<uint32_t>(owner->native.numerical_status);
  *evaluations = owner->native.evaluations;
  return COSMO_OK;
}
extern "C" uint32_t cosmo_bao_result_destroy(cosmo_bao_result *owner) {
  delete owner;
  return COSMO_OK;
}
