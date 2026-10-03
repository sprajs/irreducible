#include "irred/abi.h"
#include "irred/bao_thermal.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <new>
#include <optional>
#include <string_view>

// These declarations are schema-owned. Registration/generation is coordinated
// separately; this body introduces no parallel C layout or scientific equation.
namespace {
namespace b = irred::bao;
namespace c = irred::cosmology;
namespace n = irred::numerics;
namespace s = irred::statistics;
using Charge = irred::detail::PayloadAccounting;
constexpr size_t rows_cap = 64, models_cap = 16, species_cap = 16,
                 text_cap = 65536, bytes_cap = 64 * 1024 * 1024;
struct Slice { size_t offset = 0, length = 0; };
template <class T> bool aligned(const T *p) noexcept {
  return p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0;
}
template <class T> bool header(const T *p) noexcept {
  return aligned(p) && p->struct_size == sizeof(T) &&
         p->abi_version == IRRED_ABI_VERSION;
}
bool text(const irred_bytes &v, size_t &left) noexcept {
  if (!v.data || !v.length || v.length > 256 || v.length > left) return false;
  left -= v.length;
  return true;
}
bool ids(const irred_strings &v, size_t &left) noexcept {
  if (v.length > rows_cap || v.byte_length != v.length * sizeof(irred_bytes) ||
      (v.length && !aligned(v.data))) return false;
  for (size_t i = 0; i < v.length; ++i)
    if (!text(v.data[i], left)) return false;
  return true;
}
bool same(irred_bytes a, irred_bytes d) noexcept {
  return a.length == d.length && std::equal(a.data, a.data + a.length, d.data);
}
bool buffer(const irred_f64_buffer &v, size_t maximum) noexcept {
  return header(&v) && v.element_type == 2 && !v.reserved &&
         v.length <= maximum && v.byte_length == v.length * sizeof(double) &&
         (!v.length || aligned(v.data));
}
irred_bytes bytes(std::string_view v) noexcept {
  return {reinterpret_cast<const uint8_t *>(v.data()), v.size()};
}
irred_f64_buffer doubles(std::span<const double> v) noexcept {
  return {sizeof(irred_f64_buffer), IRRED_ABI_VERSION, 2, 0,
          v.data(), v.size(), v.size_bytes()};
}
irred_output_state state(const b::OutputState &v) noexcept {
  return {static_cast<uint32_t>(v.availability), static_cast<uint32_t>(v.status),
          static_cast<uint32_t>(v.numerical_status), 0};
}
std::string copy(irred_bytes v) {
  return {reinterpret_cast<const char *>(v.data), static_cast<size_t>(v.length)};
}
// Pre-admit growth/capacity/terminator overlap conservatively, then verify
// actual retained capacities before any scientific call. Allocator/RSS omitted.
void string_bound(Charge &q, irred_bytes v) noexcept {
  q.add(4 * static_cast<size_t>(v.length) + 64, 1);
}
void model_payload(Charge &q, const c::ThermalObservableRequest &m) noexcept {
  q.vector(m.model.species);
  q.string(m.drag_origin); q.string(m.source_origin);
}
} // namespace
struct irred_bao_thermal_result {
  // No default Gaussian metadata: both optional scientific owners disengaged.
  std::optional<b::PreparedDensity> observation;
  std::optional<b::ThermalDensityBatch> batch;
  std::vector<char> redshift_pool, model_id_pool;
  size_t redshift_length = 0;
  std::array<Slice, models_cap> model_ids{};
  std::vector<c::ThermalObservableRequest> models;
  // Destroy descriptor tables before their borrow targets (reverse order).
  std::vector<irred_bytes> ids;
  std::vector<irred_bao_query> queries;
  std::vector<irred_bao_thermal_species> species;
  std::vector<irred_bao_thermal_row> rows;
  irred_bao_thermal_view view{};
  irred_bao_thermal_result() = default;
  irred_bao_thermal_result(const irred_bao_thermal_result &) = delete;
  irred_bao_thermal_result &operator=(const irred_bao_thermal_result &) = delete;
  irred_bao_thermal_result(irred_bao_thermal_result &&) = delete;
  irred_bao_thermal_result &operator=(irred_bao_thermal_result &&) = delete;
};
namespace {
std::optional<size_t> retained(const irred_bao_thermal_result &r) noexcept {
  Charge q(sizeof(r));
  if (r.observation) q.embedded(r.observation->retained_payload_bound(), sizeof(b::PreparedDensity));
  q.vector(r.redshift_pool); q.vector(r.model_id_pool); q.vector(r.models);
  for (const auto &m : r.models) model_payload(q, m);
  q.vector(r.ids); q.vector(r.queries); q.vector(r.species); q.vector(r.rows);
  if (r.batch) {
    q.vector(r.batch->slots);
    for (const auto &x : r.batch->slots) {
      model_payload(q, x.source); q.vector(x.predictions); q.vector(x.residuals);
    }
  }
  return q.result();
}
void source_views(irred_bao_thermal_result &r) {
  if (!r.observation || r.observation->source().queries.empty()) return;
  const auto &in = r.observation->source();
  r.ids.reserve(in.ordered_ids.size()); r.queries.reserve(in.queries.size());
  for (const auto &id : in.ordered_ids) r.ids.push_back(bytes(id));
  for (const auto &x : in.queries)
    r.queries.push_back({x.z, static_cast<uint32_t>(x.observable), 0});
  auto &v = r.view;
  v.source_available = 1;
  auto &a = v.source;
  a.struct_size = sizeof(a); a.abi_version = IRRED_ABI_VERSION;
  a.role = static_cast<uint32_t>(in.role); a.covariance_unit = static_cast<uint32_t>(in.covariance_unit);
  a.queries = r.queries.data(); a.query_count = r.queries.size();
  a.query_byte_length = r.queries.size() * sizeof(irred_bao_query);
  a.observed = doubles(in.observed); a.covariance = doubles(in.covariance);
  a.ordered_ids = {r.ids.data(), r.ids.size(), r.ids.size() * sizeof(irred_bytes)};
  a.covariance_axis_ids = a.ordered_ids;
  a.table_identity = bytes(in.table_identity); a.covariance_identity = bytes(in.covariance_identity);
  a.ordering_provenance = bytes(in.ordering_provenance);
  a.calibration_provenance = bytes(in.calibration_provenance);
  a.dependence_provenance = bytes(in.dependence_provenance);
  a.redshift_convention = {reinterpret_cast<const uint8_t *>(r.redshift_pool.data()), r.redshift_length};
  const auto &meta = r.observation->metadata();
  v.actual_source_arithmetic_id = bytes(meta.arithmetic_id);
  v.source_measure = bytes(meta.measure); v.source_semantics = bytes(meta.source_semantics);
}
void row_views(irred_bao_thermal_result &r) {
  if (!r.batch) return;
  size_t count = 0;
  for (const auto &x : r.batch->slots) count += x.source.model.species.size();
  r.species.reserve(count); r.rows.reserve(r.batch->slots.size());
  for (const auto &x : r.batch->slots)
    for (const auto &sp : x.source.model.species)
      r.species.push_back({sizeof(irred_bao_thermal_species), IRRED_ABI_VERSION, 0, 0,
                          sp.mass_ev, sp.temperature_today_kelvin, sp.statistical_weight});
  size_t offset = 0;
  for (size_t i = 0; i < r.batch->slots.size(); ++i) {
    const auto &x = r.batch->slots[i]; const auto &m = x.source.model;
    irred_bao_thermal_row v{};
    v.struct_size = sizeof(v); v.abi_version = IRRED_ABI_VERSION; v.model_index = i;
    auto &a = v.source;
    a.struct_size = sizeof(a); a.abi_version = IRRED_ABI_VERSION;
    const auto id = r.model_ids[i];
    a.id = {reinterpret_cast<const uint8_t *>(r.model_id_pool.data() + id.offset), id.length};
    a.h0_km_s_mpc = m.h0_km_s_mpc;
    a.physical_baryon_density = m.physical_baryon_density;
    a.physical_cdm_density = m.physical_cdm_density; a.tcmb_kelvin = m.tcmb_kelvin;
    a.physical_massless_nonphoton_density = m.physical_massless_nonphoton_density;
    a.species_count = m.species.size(); a.species_byte_length = m.species.size() * sizeof(irred_bao_thermal_species);
    a.species = m.species.empty() ? nullptr : r.species.data() + offset;
    offset += m.species.size();
    a.z_drag = x.source.z_drag; a.drag_origin = bytes(x.source.drag_origin); a.source_origin = bytes(x.source.source_origin);
    v.preparation_status = static_cast<uint32_t>(x.preparation_status);
    v.numerical_status = static_cast<uint32_t>(x.numerical_status);
    v.prediction_state = state(x.predictions_state); v.residual_state = state(x.residuals_state);
    v.predictions = doubles(x.predictions_state.availability == c::Availability::available
                                 ? std::span<const double>(x.predictions) : std::span<const double>{});
    v.residuals = doubles(x.residuals_state.availability == c::Availability::available
                               ? std::span<const double>(x.residuals) : std::span<const double>{});
    auto &d = v.density;
    d.struct_size = sizeof(d); d.abi_version = IRRED_ABI_VERSION; d.state = state(x.density_state);
    if (x.projection_log_density_error_estimate) {
      d.has_projection_estimate = 1; d.projection_log_density_error_estimate = *x.projection_log_density_error_estimate;
    }
    if (x.result && x.density_state.availability == c::Availability::available) {
      d.quadratic = x.result->quadratic; d.log_determinant = x.result->log_determinant;
      d.log_normalization = x.result->normalization; d.log_density = x.result->density.log_value;
      d.backward_residual = x.result->backward_residual;
      d.estimated_forward_sensitivity = x.result->estimated_forward_sensitivity;
    }
    v.callbacks = x.callbacks; v.preparation_callbacks = x.preparation_callbacks;
    v.outer_callbacks = x.outer_callbacks; v.momentum_callbacks = x.momentum_callbacks;
    r.rows.push_back(v);
  }
  auto &v = r.view;
  v.rows = r.rows.data(); v.model_count = r.rows.size(); v.row_byte_length = r.rows.size() * sizeof(irred_bao_thermal_row);
  v.callbacks = r.batch->callbacks; v.preparation_callbacks = r.batch->preparation_callbacks;
  v.outer_callbacks = r.batch->outer_callbacks; v.momentum_callbacks = r.batch->momentum_callbacks;
}
} // namespace
extern "C" uint32_t irred_bao_thermal_evaluate(
    const irred_bao_thermal_source *a, const irred_bao_thermal_batch *batch,
    const irred_bao_thermal_policy *p, irred_bao_thermal_result **out) {
  if (!aligned(out)) return IRRED_INVALID_INPUT;
  *out = nullptr;
  if (!aligned(a) || !aligned(batch) || !aligned(p)) return IRRED_INVALID_INPUT;
  if (a->abi_version != IRRED_ABI_VERSION || batch->abi_version != IRRED_ABI_VERSION || p->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (!header(a) || !header(batch) || !header(p) || p->arithmetic > 1 || !p->requested || (p->requested & ~7u) ||
      p->maximum_rows > rows_cap || p->maximum_matrix_elements > rows_cap * rows_cap ||
      p->maximum_models > models_cap || p->maximum_species_per_model > species_cap ||
      p->maximum_string_bytes > text_cap || p->maximum_output_array_elements > 2048 ||
      p->maximum_native_bytes > bytes_cap || p->maximum_preparation_native_bytes > p->maximum_native_bytes ||
      p->maximum_evaluation_native_bytes > p->maximum_native_bytes || p->maximum_total_callbacks > 200000000 ||
      !std::isfinite(p->maximum_forward_sensitivity) || p->maximum_forward_sensitivity <= 0 ||
      !std::isfinite(p->maximum_projection_log_density_error) || p->maximum_projection_log_density_error <= 0)
    return IRRED_INVALID_INPUT;
  const auto &pp = p->predictions;
  if (!header(&pp) || pp.reserved || pp.momentum_method || pp.maximum_depth > 64 || pp.thermal_maximum_depth > 64 ||
      pp.maximum_callbacks_per_point > 20000000 || pp.maximum_total_callbacks > 100000000 ||
      pp.thermal_maximum_callbacks_per_evaluation > 200000 || pp.thermal_maximum_total_callbacks > 100000000 ||
      pp.maximum_native_bytes > bytes_cap || pp.thermal_maximum_native_bytes > bytes_cap)
    return IRRED_INVALID_INPUT;
  for (double x : {pp.absolute_tolerance_mpc, pp.relative_tolerance, pp.absolute_tolerance_ratio,
                   pp.relative_tolerance_ratio, pp.thermal_absolute_tolerance, pp.thermal_relative_tolerance})
    if (!std::isfinite(x) || x <= 0) return IRRED_INVALID_INPUT;
  if (a->query_count > rows_cap || batch->model_count > models_cap) return IRRED_INVALID_INPUT;
  const size_t nr = a->query_count, nm = batch->model_count;
  if (!nr || nr > rows_cap || !nm || nm > models_cap ||
      a->query_byte_length != nr * sizeof(irred_bao_query) || !aligned(a->queries) ||
      batch->model_byte_length != nm * sizeof(irred_bao_thermal_model) || !aligned(batch->models) ||
      !buffer(a->observed, rows_cap) || a->observed.length != nr ||
      !buffer(a->covariance, rows_cap * rows_cap) || a->covariance.length != nr * nr)
    return IRRED_INVALID_INPUT;
  size_t left = text_cap;
  if (!ids(a->ordered_ids, left) || !ids(a->covariance_axis_ids, left) ||
      a->ordered_ids.length != nr || a->covariance_axis_ids.length != nr) return IRRED_INVALID_INPUT;
  for (const auto *t : {&a->table_identity, &a->covariance_identity, &a->ordering_provenance,
                       &a->calibration_provenance, &a->dependence_provenance, &a->redshift_convention})
    if (!text(*t, left)) return IRRED_INVALID_INPUT;
  for (size_t i = 0; i < nr; ++i) {
    if (a->queries[i].reserved || a->queries[i].observable > 2) return IRRED_INVALID_INPUT;
    for (size_t j = 0; j < i; ++j)
      if (same(a->ordered_ids.data[i], a->ordered_ids.data[j])) return IRRED_INVALID_INPUT;
  }
  size_t ns = 0, maximum_species = 0, copied_origins = 0, id_chars = 0;
  Charge input_models(2 * nm * sizeof(c::ThermalObservableRequest));
  for (size_t i = 0; i < nm; ++i) {
    const auto &m = batch->models[i];
    if (!header(&m) || m.reserved0 || m.reserved1 || m.species_count > species_cap ||
        m.species_byte_length != m.species_count * sizeof(irred_bao_thermal_species) ||
        (m.species_count && !aligned(m.species)) || !text(m.id, left) ||
        !text(m.drag_origin, left) || !text(m.source_origin, left)) return IRRED_INVALID_INPUT;
    for (size_t j = 0; j < i; ++j) if (same(m.id, batch->models[j].id)) return IRRED_INVALID_INPUT;
    for (size_t j = 0; j < m.species_count; ++j)
      if (!header(&m.species[j]) || m.species[j].reserved0 || m.species[j].reserved1) return IRRED_INVALID_INPUT;
    ns += m.species_count; maximum_species = std::max(maximum_species, size_t(m.species_count));
    id_chars += m.id.length + 1;
    copied_origins += std::max(size_t(32), size_t(m.drag_origin.length) + 1) + std::max(size_t(32), size_t(m.source_origin.length) + 1);
    input_models.add(2 * m.species_count, sizeof(c::ThermalPhysicalSpecies));
    string_bound(input_models, m.drag_origin); string_bound(input_models, m.source_origin);
  }
  if (p->maximum_native_bytes < sizeof(irred_bao_thermal_result)) return IRRED_BAO_THERMAL_QUOTA_REFUSED;
  try {
    auto r = std::make_unique<irred_bao_thermal_result>(); auto &v = r->view;
    v.struct_size = sizeof(v); v.abi_version = IRRED_ABI_VERSION;
    v.numerical_status = static_cast<uint32_t>(n::Status::invalid_input);
    v.density_status = static_cast<uint32_t>(s::DensityStatus::invalid_input);
    v.requested = p->requested; v.arithmetic_requested = p->arithmetic;
    v.thermal_density_id = bytes(b::thermal_density_id);
    v.thermal_equation_id = bytes(c::thermal_observables_equation_id);
    v.physical_mapping_id = bytes(c::thermal_physical_mapping_id);
    v.gaussian_density_id = bytes("normalized full Gaussian ratio-coordinate density");
    v.retained_payload_bytes = v.required_global_peak_bytes = sizeof(*r);
    auto finish = [&](s::DensityStatus ds, n::Status status) {
      v.density_status = static_cast<uint32_t>(ds); v.numerical_status = static_cast<uint32_t>(status);
      const auto actual = retained(*r);
      if (!actual || *actual > p->maximum_native_bytes) {
        // Views are reset before releasing their targets. No replacement Gaussian.
        r->rows.clear(); r->species.clear(); r->ids.clear(); r->queries.clear();
        r->batch.reset(); r->observation.reset();
        std::vector<irred_bao_thermal_row>().swap(r->rows);
        std::vector<irred_bao_thermal_species>().swap(r->species);
        std::vector<irred_bytes>().swap(r->ids); std::vector<irred_bao_query>().swap(r->queries);
        std::vector<char>().swap(r->redshift_pool); std::vector<char>().swap(r->model_id_pool);
        std::vector<c::ThermalObservableRequest>().swap(r->models);
        v.source = {}; v.source_available = 0; v.rows = nullptr; v.model_count = v.row_byte_length = 0;
        v.actual_source_arithmetic_id = {}; v.source_measure = {}; v.source_semantics = {};
        v.density_status = static_cast<uint32_t>(s::DensityStatus::invalid_input);
        v.numerical_status = static_cast<uint32_t>(n::Status::work_limit);
        v.retained_payload_bytes = sizeof(*r);
      } else v.retained_payload_bytes = *actual;
      *out = r.release(); return IRRED_OK;
    };
    const size_t requested_arrays = nr * nm * (bool(p->requested & 2u) + bool(p->requested & 4u));
    if (nr > p->maximum_rows || nr * nr > p->maximum_matrix_elements || nm > p->maximum_models ||
        maximum_species > p->maximum_species_per_model || text_cap - left > p->maximum_string_bytes ||
        requested_arrays > p->maximum_output_array_elements)
      return finish(s::DensityStatus::invalid_input, n::Status::work_limit);
    if (a->role != 1 || a->covariance_unit != 0) return finish(s::DensityStatus::incompatible_metadata, n::Status::invalid_input);
    for (size_t i = 0; i < nr; ++i)
      if (!same(a->ordered_ids.data[i], a->covariance_axis_ids.data[i]))
        return finish(s::DensityStatus::incompatible_metadata, n::Status::invalid_input);
    Charge source(0), identity(0), wrapper(sizeof(*r));
    source.add(2 * nr, sizeof(b::Query) + sizeof(double) + sizeof(std::string));
    source.add(2 * nr * nr, sizeof(double));
    for (size_t i = 0; i < nr; ++i) { string_bound(source, a->ordered_ids.data[i]); string_bound(identity, a->ordered_ids.data[i]); }
    for (const auto *t : {&a->table_identity, &a->covariance_identity, &a->ordering_provenance,
                         &a->calibration_provenance, &a->dependence_provenance}) { string_bound(source, *t); string_bound(identity, *t); }
    wrapper.add(2 * (a->redshift_convention.length + 1), 1);
    wrapper.add(2 * nr, sizeof(irred_bytes) + sizeof(irred_bao_query));
    auto src = source.result(), ident = identity.result(), wrap = wrapper.result();
    auto scientific = src && ident ? b::preparation_payload_bound(nr, *src, *ident) : std::nullopt;
    Charge prep(wrap.value_or(SIZE_MAX)); if (scientific) prep.add(*scientific, 1); else prep.add(SIZE_MAX, 2);
    const auto peak = prep.result();
    if (!peak) return finish(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.preparation_peak_bound_bytes = *peak; v.required_global_peak_bytes = *peak;
    if (*peak > p->maximum_native_bytes || *peak > p->maximum_preparation_native_bytes)
      return finish(s::DensityStatus::invalid_input, n::Status::work_limit);
    r->redshift_length = a->redshift_convention.length;
    r->redshift_pool.resize(r->redshift_length + 1);
    std::copy_n(a->redshift_convention.data, r->redshift_length, r->redshift_pool.data());
    r->ids.reserve(nr); r->queries.reserve(nr);
    b::DensityInput in;
    in.queries.reserve(nr); in.observed.assign(a->observed.data, a->observed.data + nr);
    in.covariance.assign(a->covariance.data, a->covariance.data + nr * nr);
    in.ordered_ids.reserve(nr);
    for (size_t i = 0; i < nr; ++i) {
      in.queries.push_back({a->queries[i].z, static_cast<b::Observable>(a->queries[i].observable)});
      in.ordered_ids.push_back(copy(a->ordered_ids.data[i]));
    }
    in.role = b::RowRole::synthetic_control; in.covariance_unit = b::CovarianceUnit::dimensionless_ratio_squared;
    in.table_identity = copy(a->table_identity); in.covariance_identity = copy(a->covariance_identity);
    in.ordering_provenance = copy(a->ordering_provenance); in.calibration_provenance = copy(a->calibration_provenance);
    in.dependence_provenance = copy(a->dependence_provenance);
    const auto actual_source = b::retained_source_payload_bound(in);
    if (!actual_source || *actual_source > *src || r->redshift_pool.capacity() > 2 * (r->redshift_length + 1) ||
        r->ids.capacity() > 2 * nr || r->queries.capacity() > 2 * nr)
      return finish(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.phase = 1; v.source_prepare_call_attempted = 1;
    r->observation.emplace(b::prepare_density(std::move(in), {size_t(p->maximum_rows), size_t(p->maximum_matrix_elements),
      size_t(p->maximum_string_bytes), size_t(p->maximum_preparation_native_bytes - *wrap), p->maximum_forward_sensitivity,
      static_cast<n::Arithmetic>(p->arithmetic)}));
    source_views(*r);
    v.source_factor_completed = r->observation->status() == s::DensityStatus::finite;
    if (!v.source_factor_completed) return finish(r->observation->status(), r->observation->numerical_status());
    v.phase = 2;
    const auto prefix = retained(*r), model_bound = input_models.result();
    auto native = b::thermal_density_payload_bound(nm, nr, copied_origins, ns, maximum_species, p->requested);
    Charge eval(prefix.value_or(SIZE_MAX));
    if (!model_bound || !native) eval.add(SIZE_MAX, 2);
    else { eval.add(*model_bound, 1); eval.add(*native, 1); }
    eval.add(2 * id_chars, 1); eval.add(2 * nm, sizeof(irred_bao_thermal_row));
    eval.add(2 * ns, sizeof(irred_bao_thermal_species));
    const auto ev = eval.result();
    if (!ev) return finish(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.evaluation_peak_bound_bytes = *ev; v.required_global_peak_bytes = std::max(*peak, *ev);
    if (*ev > p->maximum_native_bytes || *ev > p->maximum_evaluation_native_bytes)
      return finish(s::DensityStatus::invalid_input, n::Status::work_limit);
    r->model_id_pool.resize(id_chars); r->models.reserve(nm); r->rows.reserve(nm); r->species.reserve(ns);
    size_t offset = 0;
    for (size_t i = 0; i < nm; ++i) {
      const auto &m = batch->models[i]; r->model_ids[i] = {offset, size_t(m.id.length)};
      std::copy_n(m.id.data, m.id.length, r->model_id_pool.data() + offset); offset += m.id.length + 1;
      c::ThermalObservableRequest req{{m.h0_km_s_mpc, m.physical_baryon_density, m.physical_cdm_density,
        m.tcmb_kelvin, m.physical_massless_nonphoton_density, {}}, m.z_drag, copy(m.drag_origin), copy(m.source_origin)};
      req.model.species.reserve(m.species_count);
      for (size_t j = 0; j < m.species_count; ++j) {
        const auto &sp = m.species[j]; req.model.species.push_back({sp.mass_ev, sp.temperature_today_kelvin, sp.statistical_weight});
      }
      r->models.push_back(std::move(req));
    }
    Charge actual_models(0); actual_models.vector(r->models);
    for (const auto &m : r->models) model_payload(actual_models, m);
    const auto allocated_models = actual_models.result();
    if (!allocated_models || *allocated_models > *model_bound || r->model_id_pool.capacity() > 2 * id_chars ||
        r->rows.capacity() > 2 * nm || r->species.capacity() > 2 * ns)
      return finish(s::DensityStatus::invalid_input, n::Status::work_limit);
    b::ThermalDensityPolicy policy;
    policy.maximum_models = p->maximum_models; policy.maximum_queries = p->maximum_rows;
    policy.maximum_string_bytes = p->maximum_string_bytes;
    // The kernel helper excludes retained source/borrowed owned input. Charge
    // those wrapper owners before assigning the original remaining subcap.
    const size_t wrapper_eval = *ev - *native;
    policy.maximum_native_bytes = std::min(size_t(p->maximum_native_bytes), size_t(p->maximum_evaluation_native_bytes)) - wrapper_eval;
    policy.maximum_total_callbacks = p->maximum_total_callbacks; policy.maximum_forward_sensitivity = p->maximum_forward_sensitivity;
    policy.maximum_projection_log_density_error = p->maximum_projection_log_density_error;
    policy.arithmetic = static_cast<n::Arithmetic>(p->arithmetic); policy.requested = p->requested;
    auto &pr = policy.predictions;
    pr.absolute_tolerance_mpc = pp.absolute_tolerance_mpc; pr.relative_tolerance = pp.relative_tolerance;
    pr.absolute_tolerance_ratio = pp.absolute_tolerance_ratio; pr.relative_tolerance_ratio = pp.relative_tolerance_ratio;
    pr.maximum_callbacks_per_point = pp.maximum_callbacks_per_point; pr.maximum_total_callbacks = pp.maximum_total_callbacks;
    pr.maximum_depth = pp.maximum_depth; pr.maximum_points = nr;
    pr.maximum_native_bytes = std::min(size_t(pp.maximum_native_bytes), policy.maximum_native_bytes);
    pr.thermal.absolute_tolerance = pp.thermal_absolute_tolerance; pr.thermal.relative_tolerance = pp.thermal_relative_tolerance;
    pr.thermal.maximum_callbacks_per_evaluation = pp.thermal_maximum_callbacks_per_evaluation;
    pr.thermal.maximum_total_callbacks = pp.thermal_maximum_total_callbacks; pr.thermal.maximum_depth = pp.thermal_maximum_depth;
    pr.thermal.maximum_points = nr; pr.thermal.maximum_species = p->maximum_species_per_model;
    pr.thermal.maximum_native_bytes = std::min(size_t(pp.thermal_maximum_native_bytes), policy.maximum_native_bytes);
    pr.thermal.momentum_method = c::ThermalMomentumMethod::direct_adaptive;
    v.phase = 2; v.thermal_batch_call_attempted = 1;
    r->batch.emplace(r->observation->evaluate_thermal(r->models, policy));
    row_views(*r); v.phase = 3;
    return finish(r->batch->status, r->batch->numerical_status);
  } catch (const std::bad_alloc &) { return IRRED_ALLOCATION_FAILURE; }
    catch (...) { return IRRED_EXCEPTION; }
}
extern "C" uint32_t irred_bao_thermal_result_view(const irred_bao_thermal_result *r, irred_bao_thermal_view *out) {
  if (!aligned(r) || !aligned(out)) return IRRED_INVALID_INPUT;
  *out = r->view; return IRRED_OK;
}
extern "C" uint32_t irred_bao_thermal_result_destroy(irred_bao_thermal_result *r) {
  if (r && !aligned(r)) return IRRED_INVALID_INPUT;
  delete r; return IRRED_OK;
}
