#include "irred/abi.h"
#include "irred/bao.hpp"
#include "payload_accounting.hpp"
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>
namespace {
namespace bao = irred::bao;
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
irred_bytes text(std::string_view s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
irred_strings strings(const std::vector<irred_bytes> &v) {
  return {v.data(), v.size(), v.size() * sizeof(irred_bytes)};
}
irred_f64_buffer doubles(const std::vector<double> &v) {
  return {sizeof(irred_f64_buffer), IRRED_ABI_VERSION, 2, 0, v.data(), v.size(),
          v.size() * sizeof(double)};
}
irred_u64_buffer indices(const std::vector<size_t> &v) {
  static_assert(sizeof(size_t) == sizeof(uint64_t));
  return {reinterpret_cast<const uint64_t *>(v.data()), v.size(),
          v.size() * sizeof(uint64_t)};
}
irred_output_state state(const bao::OutputState &s) {
  return {(uint32_t)s.availability, (uint32_t)s.status,
          (uint32_t)s.numerical_status, 0};
}
bool precision(uint32_t a, double b) {
  return a <= 1 && std::isfinite(b) && b > 0;
}
bool projection(const irred_projection_policy &p) {
  if (p.struct_size != sizeof(p) || p.abi_version != IRRED_ABI_VERSION ||
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
c::EvaluationPolicy native_projection(const irred_projection_policy &p,
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
bool model_descriptor(const irred_bao_model &m) {
  const auto &e = m.expansion;
  const auto &b = e.parameters;
  const size_t n = e.model == 2 ? 3 : e.model == 3 ? 5 : 1;
  return m.struct_size == sizeof(m) && m.abi_version == IRRED_ABI_VERSION &&
         !m.geometry && !m.reserved && e.struct_size == sizeof(e) &&
         e.abi_version == IRRED_ABI_VERSION && !e.reserved && e.model <= 3 &&
         b.struct_size == sizeof(b) && b.abi_version == IRRED_ABI_VERSION &&
         b.element_type == 2 && !b.reserved && b.length == n &&
         descriptor(b.data, b.length, b.byte_length);
}
c::ExpansionSpec expansion(uint32_t model, const double *p) {
  switch (model) {
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
bool buffer(const irred_f64_buffer &b) {
  return b.struct_size == sizeof(b) && b.abi_version == IRRED_ABI_VERSION &&
         b.element_type == 2 && !b.reserved &&
         descriptor(b.data, b.length, b.byte_length);
}
bool bytes(const irred_bytes &b) { return !b.length || b.data; }
std::string_view value(const irred_bytes &b) {
  return {reinterpret_cast<const char *>(b.data), b.length};
}
struct Core {
  bao::PreparedDensity native;
  uint32_t arithmetic = 0;
  std::string redshift, ruler;
  std::vector<irred_bytes> ids;
  std::vector<irred_bao_query> queries;
  irred_bao_view view{};
  void initialize() {
    const auto &s = native.source();
    ids.reserve(s.ordered_ids.size());
    queries.reserve(s.queries.size());
    for (const auto &id : s.ordered_ids)
      ids.push_back(text(id));
    for (const auto &q : s.queries)
      queries.push_back({q.z, (uint32_t)q.observable, 0});
    view.struct_size = sizeof(view);
    view.abi_version = IRRED_ABI_VERSION;
    view.status = (uint32_t)native.status();
    view.numerical_status = (uint32_t)native.numerical_status();
    view.arithmetic = arithmetic;
    auto &v = view.source;
    v.struct_size = sizeof(v);
    v.abi_version = IRRED_ABI_VERSION;
    v.role = (uint32_t)s.role;
    v.covariance_unit = (uint32_t)s.covariance_unit;
    v.queries = queries.data();
    v.query_count = queries.size();
    v.query_byte_length = queries.size() * sizeof(irred_bao_query);
    v.observed = doubles(s.observed);
    v.covariance = doubles(s.covariance);
    v.ordered_ids = strings(ids);
    v.covariance_axis_ids = strings(ids);
    v.table_identity = text(s.table_identity);
    v.covariance_identity = text(s.covariance_identity);
    v.ordering_provenance = text(s.ordering_provenance);
    v.calibration_provenance = text(s.calibration_provenance);
    v.dependence_provenance = text(s.dependence_provenance);
    v.redshift_convention = text(redshift);
    v.ruler_convention = text(ruler);
    view.arithmetic_id =
        text(arithmetic ? "F02/longdouble-cpu/v1" : "F02/binary64-legacy/v1");
    view.density_id = text("F03/normalized-Gaussian/v1");
    view.equation_id = text(bao::equation_id);
    view.ruler_convention_id = text(bao::ruler_convention_id);
    irred::detail::PayloadAccounting b(sizeof(Core));
    b.embedded(native.retained_payload_bound(), sizeof(native));
    b.vector(ids);
    b.vector(queries);
    b.string(redshift);
    b.string(ruler);
    auto total = b.result();
    if (!total)
      throw std::bad_alloc();
    view.retained_bytes = *total;
  }
};
} // namespace
struct irred_bao {
  std::shared_ptr<Core> core;
};
struct irred_bao_result {
  std::shared_ptr<Core> core;
  bao::DensityBatch native;
  std::vector<std::vector<double>> parameters;
  std::vector<irred_bao_row> rows;
};
extern "C" uint32_t
irred_bao_prepare(const irred_bao_source *s,
                          const irred_bao_preparation_policy *p,
                          irred_bao **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(s) || !aligned(p))
    return IRRED_INVALID_INPUT;
  if (s->abi_version != IRRED_ABI_VERSION ||
      p->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (s->struct_size != sizeof(*s) || p->struct_size != sizeof(*p) ||
      p->reserved ||
      !precision(p->arithmetic, p->maximum_forward_sensitivity) ||
      p->maximum_queries > 4096 || p->maximum_matrix_elements > SIZE_MAX ||
      p->maximum_string_bytes > SIZE_MAX ||
      p->maximum_native_bytes > SIZE_MAX || s->query_count > 4096 ||
      !descriptor(s->queries, s->query_count, s->query_byte_length) ||
      !buffer(s->observed) || !buffer(s->covariance) ||
      !descriptor(s->ordered_ids.data, s->ordered_ids.length,
                  s->ordered_ids.byte_length) ||
      !descriptor(s->covariance_axis_ids.data, s->covariance_axis_ids.length,
                  s->covariance_axis_ids.byte_length))
    return IRRED_INVALID_INPUT;
  const size_t n = s->query_count;
  if (!n || n > p->maximum_queries || n > p->maximum_matrix_elements / n ||
      s->observed.length != n || s->covariance.length != n * n ||
      s->ordered_ids.length != n || s->covariance_axis_ids.length != n ||
      (s->role > 1) || s->covariance_unit)
    return IRRED_INVALID_INPUT;
  try {
    irred::detail::PayloadAccounting source(0), identity(0),
        wrapper(sizeof(Core) + sizeof(irred_bao));
    source.add(n, sizeof(bao::Query) + sizeof(double) + sizeof(std::string));
    source.add(n * n, sizeof(double));
    wrapper.add(n, sizeof(irred_bytes) + sizeof(irred_bao_query));
    size_t chars = 0;
    auto string_payload = [](const irred_bytes &b) -> std::optional<size_t> {
      if (!bytes(b) || b.length == SIZE_MAX)
        return {};
      return std::max<size_t>(b.length, std::string{}.capacity()) + 1;
    };
    for (size_t i = 0; i < n; ++i) {
      const auto &id = s->ordered_ids.data[i],
                 &axis = s->covariance_axis_ids.data[i];
      auto size = string_payload(id);
      if (!size || !bytes(axis) || !id.length || axis.length != id.length ||
          !add(chars, id.length, 1) || chars > p->maximum_string_bytes)
        return IRRED_INVALID_INPUT;
      // Bound logical text before inspecting declared ID/axis contents.
      if (value(id) != value(axis))
        return IRRED_INVALID_INPUT;
      source.add(*size, 1);
      identity.add(*size, 1);
      if (s->queries[i].reserved || s->queries[i].observable > 2)
        return IRRED_INVALID_INPUT;
    }
    for (const auto *v :
         {&s->table_identity, &s->covariance_identity, &s->ordering_provenance,
          &s->calibration_provenance, &s->dependence_provenance}) {
      auto size = string_payload(*v);
      if (!size || !v->length || !add(chars, v->length, 1) ||
          chars > p->maximum_string_bytes)
        return IRRED_INVALID_INPUT;
      source.add(*size, 1);
      identity.add(*size, 1);
    }
    for (const auto *v : {&s->redshift_convention, &s->ruler_convention}) {
      auto size = string_payload(*v);
      if (!size || !v->length || !add(chars, v->length, 1) ||
          chars > p->maximum_string_bytes)
        return IRRED_INVALID_INPUT;
      wrapper.add(*size, 1);
    }
    auto sb = source.result(), ib = identity.result(), wb = wrapper.result();
    auto peak =
        sb && ib ? bao::preparation_payload_bound(n, *sb, *ib) : std::nullopt;
    if (chars > p->maximum_string_bytes || !peak || !wb ||
        *wb > p->maximum_native_bytes || *peak > p->maximum_native_bytes - *wb)
      return IRRED_INVALID_INPUT;
    bao::DensityInput input;
    input.queries.reserve(n);
    input.ordered_ids.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      input.queries.push_back(
          {s->queries[i].z, (bao::Observable)s->queries[i].observable});
      input.ordered_ids.emplace_back(value(s->ordered_ids.data[i]));
    }
    input.observed.assign(s->observed.data, s->observed.data + n);
    input.covariance.assign(s->covariance.data, s->covariance.data + n * n);
    input.role = (bao::RowRole)s->role;
    input.covariance_unit = (bao::CovarianceUnit)s->covariance_unit;
    input.table_identity = value(s->table_identity);
    input.covariance_identity = value(s->covariance_identity);
    input.ordering_provenance = value(s->ordering_provenance);
    input.calibration_provenance = value(s->calibration_provenance);
    input.dependence_provenance = value(s->dependence_provenance);
    auto core = std::make_shared<Core>();
    core->arithmetic = p->arithmetic;
    core->redshift = value(s->redshift_convention);
    core->ruler = value(s->ruler_convention);
    bao::PreparationPolicy policy{(size_t)p->maximum_queries,
                                  (size_t)p->maximum_matrix_elements,
                                  (size_t)p->maximum_string_bytes,
                                  (size_t)(p->maximum_native_bytes - *wb),
                                  p->maximum_forward_sensitivity,
                                  (num::Arithmetic)p->arithmetic};
    core->native = bao::prepare_density(std::move(input), policy);
    if (core->native.numerical_status() == num::Status::work_limit ||
        core->native.source().queries.empty())
      return IRRED_INVALID_INPUT;
    core->initialize();
    auto owner = std::make_unique<irred_bao>();
    owner->core = std::move(core);
    *out = owner.release();
    return IRRED_OK;
  } catch (const std::bad_alloc &) {
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t irred_bao_source_view(const irred_bao *s,
                                                  irred_bao_view *v) {
  if (aligned(v))
    *v = {};
  if (!aligned(s) || !aligned(v))
    return IRRED_INVALID_INPUT;
  *v = s->core->view;
  return IRRED_OK;
}
extern "C" uint32_t irred_bao_destroy(irred_bao *s) {
  if (!s)
    return IRRED_OK;
  if (!aligned(s))
    return IRRED_INVALID_INPUT;
  delete s;
  return IRRED_OK;
}
extern "C" uint32_t
irred_bao_evaluate(const irred_bao *s,
                           const irred_bao_batch *b,
                           const irred_bao_evaluation_policy *p,
                           irred_bao_result **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(s) || !aligned(b) || !aligned(p))
    return IRRED_INVALID_INPUT;
  if (b->abi_version != IRRED_ABI_VERSION ||
      p->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (b->struct_size != sizeof(*b) || p->struct_size != sizeof(*p) ||
      !precision(p->arithmetic, p->maximum_forward_sensitivity) ||
      !p->requested || (p->requested & ~7u) || !projection(p->projection) ||
      p->maximum_models > 65536 || p->maximum_array_elements > 1048576 ||
      p->maximum_native_bytes > SIZE_MAX || b->model_count > 65536 ||
      !descriptor(b->models, b->model_count, b->model_byte_length))
    return IRRED_INVALID_INPUT;
  const auto selected_rows = s->core->native.source().queries.size();
  if (selected_rows && b->model_count > 1048576 / selected_rows)
    return IRRED_INVALID_INPUT;
  try {
    auto result = std::make_unique<irred_bao_result>();
    if (s->core->native.status() != irred::statistics::DensityStatus::finite) {
      result->native.status = s->core->native.status();
      result->native.numerical_status = s->core->native.numerical_status();
      result->core = s->core;
      *out = result.release();
      return IRRED_OK;
    }
    if (p->arithmetic != s->core->arithmetic) {
      result->native.status =
          irred::statistics::DensityStatus::incompatible_metadata;
      result->core = s->core;
      *out = result.release();
      return IRRED_OK;
    }
    const size_t count = b->model_count,
                 n = s->core->native.source().queries.size();
    const unsigned arrays = bool(p->requested & 2u) + bool(p->requested & 4u);
    if (n && count > SIZE_MAX / n)
      return IRRED_INVALID_INPUT;
    const size_t elements = n * count;
    if (arrays && elements > SIZE_MAX / arrays)
      return IRRED_INVALID_INPUT;
    size_t used = sizeof(irred_bao_result);
    if (!add(used, count,
             sizeof(irred_bao_row) + sizeof(bao::ModelPoint) +
                 sizeof(std::vector<double>) + 5 * sizeof(double)))
      return IRRED_INVALID_INPUT;
    if (count > p->maximum_models ||
        arrays * elements > p->maximum_array_elements ||
        used > p->maximum_native_bytes) {
      result->native.status =
          irred::statistics::DensityStatus::numerical_failure;
      result->native.numerical_status = num::Status::work_limit;
      *out = result.release();
      return IRRED_OK;
    }
    for (size_t i = 0; i < count; ++i)
      if (!model_descriptor(b->models[i]))
        return IRRED_INVALID_INPUT;
    std::vector<bao::ModelPoint> models;
    models.reserve(count);
    result->parameters.reserve(count);
    result->rows.reserve(count);
    for (size_t i = 0; i < count; ++i) {
      const auto &m = b->models[i];
      result->parameters.emplace_back(m.expansion.parameters.data,
                                      m.expansion.parameters.data +
                                          m.expansion.parameters.length);
      models.emplace_back(
          expansion(m.expansion.model, m.expansion.parameters.data),
          c::FlatFLRW{}, bao::Ruler(m.h0_rd_km_s));
    }
    bao::DensityPolicy policy;
    policy.arithmetic = (num::Arithmetic)p->arithmetic;
    policy.maximum_forward_sensitivity = p->maximum_forward_sensitivity;
    policy.maximum_models = p->maximum_models;
    policy.maximum_native_bytes = p->maximum_native_bytes - used;
    policy.requested = p->requested;
    policy.observables.maximum_queries = p->projection.maximum_queries;
    policy.observables.maximum_native_bytes = policy.maximum_native_bytes;
    policy.observables.background =
        native_projection(p->projection, policy.maximum_native_bytes);
    result->native = s->core->native.evaluate(models, policy);
    result->core = s->core;
    for (size_t i = 0; i < result->native.slots.size(); ++i) {
      const auto &v = result->native.slots[i];
      irred_bao_row r{};
      r.struct_size = sizeof(r);
      r.abi_version = IRRED_ABI_VERSION;
      r.model_index = i;
      r.source = b->models[i];
      r.source.expansion.parameters.data = result->parameters[i].data();
      r.background_status = (uint32_t)v.background_status;
      r.numerical_status = (uint32_t)v.numerical_status;
      r.predictions_state = state(v.predictions_state);
      r.residuals_state = state(v.residuals_state);
      r.density.state = state(v.density_state);
      if (v.result && v.result->density.status ==
                          irred::statistics::DensityStatus::finite) {
        const auto &g = *v.result;
        r.density.quadratic = g.quadratic;
        r.density.log_determinant = g.log_determinant;
        r.density.log_normalization = g.normalization;
        r.density.log_density = g.density.log_value;
        r.density.backward_residual = g.backward_residual;
        r.density.estimated_forward_sensitivity =
            g.estimated_forward_sensitivity;
      }
      r.predictions = doubles(v.predictions);
      r.residuals = doubles(v.residuals);
      r.background_node_indices = indices(v.node_indices);
      r.callbacks = v.work.callbacks;
      r.segment_visits = v.work.segment_visits;
      r.model_id = text(v.model_id);
      r.arithmetic_id = s->core->view.arithmetic_id;
      result->rows.push_back(r);
    }
    *out = result.release();
    return IRRED_OK;
  } catch (const std::bad_alloc &) {
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t irred_bao_result_view(
    const irred_bao_result *r, const irred_bao_row **rows,
    uint64_t *count, uint32_t *status, uint32_t *numerical, uint64_t *callbacks,
    uint64_t *segments) {
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
    return IRRED_INVALID_INPUT;
  *rows = r->rows.data();
  *count = r->rows.size();
  *status = (uint32_t)r->native.status;
  *numerical = (uint32_t)r->native.numerical_status;
  *callbacks = r->native.work.callbacks;
  *segments = r->native.work.segment_visits;
  return IRRED_OK;
}
extern "C" uint32_t
irred_bao_result_source_view(const irred_bao_result *r,
                                     irred_bao_view *v,
                                     uint32_t *available) {
  if (aligned(v))
    *v = {};
  if (aligned(available))
    *available = 0;
  if (!aligned(r) || !aligned(v) || !aligned(available))
    return IRRED_INVALID_INPUT;
  if (r->core) {
    *v = r->core->view;
    *available = 1;
  }
  return IRRED_OK;
}
extern "C" uint32_t
irred_bao_result_destroy(irred_bao_result *r) {
  if (!r)
    return IRRED_OK;
  if (!aligned(r))
    return IRRED_INVALID_INPUT;
  delete r;
  return IRRED_OK;
}
