#include "payload_accounting.hpp"
#include "irred/abi.h"
#include "irred/background.hpp"
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>
namespace {
namespace c = irred::cosmology;
namespace n = irred::numerics;
template <class T> bool aligned(const T *p) {
  return p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0;
}
template <class T> bool descriptor(const T *p, uint64_t count, uint64_t bytes) {
  return count <= SIZE_MAX / sizeof(T) && bytes == count * sizeof(T) &&
         (!count || aligned(p));
}
bool add(size_t &bytes, size_t count, size_t width) {
  return irred::detail::checked_payload_add(bytes, count, width);
}
irred_bytes view(std::string_view s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
uint32_t numerical(c::Status s) {
  if (s == c::Status::ok)
    return (uint32_t)n::Status::ok;
  if (s == c::Status::work_limit)
    return (uint32_t)n::Status::work_limit;
  if (s == c::Status::unsupported_domain)
    return (uint32_t)n::Status::outside_domain;
  return (uint32_t)n::Status::invalid_input;
}
template <class T> irred_output_state state(const c::Outcome<T> &o) {
  return {(uint32_t)o.availability, (uint32_t)o.status,
          (uint32_t)o.numerical_status, 0};
}
irred_scalar_outcome scalar(const c::Outcome<double> &o) {
  return {state(o), o.value.value_or(0)};
}
c::ExpansionSpec spec(uint32_t model, const std::vector<double> &p) {
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
irred_expansion_row row(const c::Slot &s, size_t mi, size_t qi, size_t offset) {
  irred_expansion_row r{};
  r.struct_size = sizeof(r);
  r.abi_version = IRRED_ABI_VERSION;
  r.model_index = mi;
  r.query_index = qi;
  r.admission_status = (uint32_t)s.admission_status;
  if (s.node_index) {
    r.has_node = 1;
    r.node_index = offset + *s.node_index;
  }
  r.expansion.state = state(s.expansion);
  if (s.expansion.value) {
    r.expansion.E = s.expansion.value->expansion_E;
    r.expansion.H_km_s_mpc = scalar(s.expansion.value->h_km_s_mpc);
  }
  r.radial.state = state(s.radial);
  if (s.radial.value) {
    r.radial.E = s.radial.value->expansion_E;
    r.radial.integral = s.radial.value->integral;
    r.radial.error_estimate = s.radial.value->error_estimate;
  }
  r.luminosity_shape = scalar(s.luminosity_shape);
  r.clock.state = state(s.clock);
  if (s.clock.value) {
    r.clock.integral = s.clock.value->integral;
    r.clock.error_estimate = s.clock.value->error_estimate;
    r.clock.lookback_seconds = scalar(s.clock.value->lookback_seconds);
  }
  r.physical.state = state(s.physical);
  if (s.physical.value) {
    const auto &v = *s.physical.value;
    r.physical.radial_mpc = v.radial_mpc;
    r.physical.transverse_mpc = v.transverse_mpc;
    r.physical.angular_diameter_mpc = v.angular_diameter_mpc;
    r.physical.luminosity_mpc = v.luminosity_mpc;
    r.physical.volume_mpc3_per_sr_per_redshift =
        v.volume_mpc3_per_sr_per_redshift;
  }
  r.kinematics.state = state(s.kinematics);
  if (s.kinematics.value) {
    const auto &v = *s.kinematics.value;
    r.kinematics.q = v.q;
    r.kinematics.has_jerk = v.jerk.has_value();
    r.kinematics.jerk = v.jerk.value_or(0);
    r.kinematics.has_q0 = v.q0_within_piecewise_model.has_value();
    r.kinematics.q0_within_piecewise_model =
        v.q0_within_piecewise_model.value_or(0);
    r.kinematics.q_convention = (uint32_t)v.q_convention;
    r.kinematics.jerk_availability = (uint32_t)v.jerk_availability;
    r.kinematics.bin = v.bin;
  }
  return r;
}
} // namespace
struct irred_expansion_result {
  std::vector<std::vector<double>> parameters;
  std::vector<irred_expansion_model_view> models;
  std::vector<irred_expansion_request> queries;
  std::vector<irred_expansion_row> rows;
  std::vector<irred_expansion_node> nodes;
  uint32_t status = (uint32_t)c::Status::ok,
           numerical_status = (uint32_t)n::Status::ok;
  uint64_t callbacks = 0, segments = 0;
  bool source_available = false;
};
extern "C" uint32_t irred_expansion_evaluate(const irred_expansion_batch *b,
                                             const irred_expansion_policy *p,
                                             irred_expansion_result **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(b) || !aligned(p))
    return IRRED_INVALID_INPUT;
  if (b->abi_version != IRRED_ABI_VERSION ||
      p->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (b->struct_size != sizeof(*b) || p->struct_size != sizeof(*p) ||
      b->reserved || b->geometry || p->has_integration > 1 ||
      p->maximum_models > 65536 || p->maximum_queries > 65536 ||
      p->maximum_slots > 1048576 || b->model_count > 65536 ||
      b->query_count > 65536 ||
      !descriptor(b->models, b->model_count, b->model_byte_length) ||
      !descriptor(b->queries, b->query_count, b->query_byte_length))
    return IRRED_INVALID_INPUT;
  if (p->maximum_callbacks > SIZE_MAX || p->maximum_segment_visits > SIZE_MAX ||
      p->maximum_native_bytes > SIZE_MAX ||
      p->integration_max_evaluations > SIZE_MAX)
    return IRRED_INVALID_INPUT;
  if (!p->has_integration) {
    if (p->max_depth || p->absolute_tolerance != 0 ||
        p->relative_tolerance != 0 || p->integration_max_evaluations)
      return IRRED_INVALID_INPUT;
  } else if (p->max_depth > 60 || !std::isfinite(p->absolute_tolerance) ||
             !std::isfinite(p->relative_tolerance) ||
             p->absolute_tolerance < 0 || p->relative_tolerance < 0 ||
             (p->absolute_tolerance == 0 && p->relative_tolerance == 0))
    return IRRED_INVALID_INPUT;
  if (b->query_count && b->model_count > SIZE_MAX / b->query_count)
    return IRRED_INVALID_INPUT;
  const size_t total = b->model_count * b->query_count;
  bool capped = b->model_count > p->maximum_models ||
                b->query_count > p->maximum_queries || total > p->maximum_slots;
  size_t bytes = sizeof(irred_expansion_result);
  const auto workspace = c::Expansion::workspace_payload_bound(b->query_count);
  if (!workspace ||
      !add(bytes, b->model_count,
           sizeof(std::vector<double>) + 5 * sizeof(double) +
               sizeof(irred_expansion_model_view)) ||
      !add(bytes, b->query_count,
           sizeof(irred_expansion_request) + sizeof(c::Request)) ||
      !add(bytes, total,
           sizeof(irred_expansion_row) + sizeof(irred_expansion_node)) ||
      !add(bytes, 1, *workspace) || !add(bytes, 1, sizeof(c::Expansion)))
    return IRRED_INVALID_INPUT;
  capped |= bytes > p->maximum_native_bytes;
  try {
    auto owner = std::make_unique<irred_expansion_result>();
    if (capped) {
      owner->status = (uint32_t)c::Status::work_limit;
      owner->numerical_status = (uint32_t)n::Status::work_limit;
      *out = owner.release();
      return IRRED_OK;
    }
    // Only admitted descriptors are dereferenced. Model-domain errors are
    // native scientific outcomes; malformed tags/storage remain boundary
    // failures.
    for (size_t i = 0; i < b->model_count; ++i) {
      const auto &m = b->models[i];
      const auto &a = m.parameters;
      size_t count = m.model == 2 ? 3 : m.model == 3 ? 5 : 1;
      if (m.struct_size != sizeof(m) || m.abi_version != IRRED_ABI_VERSION ||
          m.reserved || m.model > 3 || a.struct_size != sizeof(a) ||
          a.abi_version != IRRED_ABI_VERSION || a.element_type != 2 ||
          a.reserved || a.length != count ||
          !descriptor(a.data, a.length, a.byte_length))
        return IRRED_INVALID_INPUT;
    }
    for (size_t i = 0; i < b->query_count; ++i) {
      const auto &q = b->queries[i];
      if (q.struct_size != sizeof(q) || q.abi_version != IRRED_ABI_VERSION ||
          q.reserved || !q.requested || (q.requested & ~63u) ||
          (q.presence_flags & ~3u) ||
          (!(q.presence_flags & 1) &&
           (q.observer_redshift != 0 || q.observer_convention)) ||
          ((q.presence_flags & 1) && q.observer_convention > 1) ||
          (!(q.presence_flags & 2) && q.h0_km_s_mpc != 0))
        return IRRED_INVALID_INPUT;
    }
    owner->parameters.reserve(b->model_count);
    owner->models.reserve(b->model_count);
    owner->queries.reserve(b->query_count);
    owner->rows.reserve(total);
    owner->nodes.reserve(total);
    std::vector<c::Request> requests;
    requests.reserve(b->query_count);
    for (size_t i = 0; i < b->query_count; ++i) {
      const auto &q = b->queries[i];
      owner->queries.push_back(q);
      std::optional<c::Observer> o;
      std::optional<c::PhysicalScale> h;
      if (q.presence_flags & 1)
        o = c::Observer{q.observer_redshift,
                        (c::Convention)q.observer_convention};
      if (q.presence_flags & 2)
        h = c::PhysicalScale(q.h0_km_s_mpc);
      requests.emplace_back(q.z_expansion, q.requested, o, h);
    }
    c::EvaluationPolicy policy;
    policy.maximum_queries = p->maximum_queries;
    policy.maximum_native_bytes = *workspace;
    if (p->has_integration)
      policy.integration = n::IntegrationPolicy{
          p->absolute_tolerance, p->relative_tolerance,
          (size_t)p->integration_max_evaluations, p->max_depth};
    for (size_t i = 0; i < b->model_count; ++i) {
      const auto &m = b->models[i];
      owner->parameters.emplace_back(m.parameters.data,
                                     m.parameters.data + m.parameters.length);
      auto expansion =
          c::prepare(spec(m.model, owner->parameters.back()), c::FlatFLRW{});
      policy.maximum_callbacks = p->maximum_callbacks - owner->callbacks;
      policy.maximum_segment_visits =
          p->maximum_segment_visits - owner->segments;
      auto native = expansion.evaluate(requests, policy);
      irred_expansion_model_view model{};
      model.struct_size = sizeof(model);
      model.abi_version = IRRED_ABI_VERSION;
      model.model_index = i;
      model.preparation_status = (uint32_t)expansion.status();
      model.evaluation_status = (uint32_t)native.status;
      model.numerical_status = numerical(native.status);
      model.row_offset = owner->rows.size();
      model.row_count = native.slots.size();
      model.source = m;
      model.source.parameters.data = owner->parameters.back().data();
      model.model_id = view(expansion.model_id());
      model.constants_id = view(c::Expansion::constants_id);
      model.radial_equation_id = view(c::Expansion::radial_equation_id);
      const auto nodeoffset = owner->nodes.size();
      for (size_t j = 0; j < native.slots.size(); ++j)
        owner->rows.push_back(row(native.slots[j], i, j, nodeoffset));
      for (const auto &node : native.nodes)
        owner->nodes.push_back({sizeof(irred_expansion_node), IRRED_ABI_VERSION,
                                i, node.z_expansion, node.work.callbacks,
                                node.work.segment_visits});
      owner->callbacks += native.work.callbacks;
      owner->segments += native.work.segment_visits;
      owner->models.push_back(model);
    }
    owner->source_available = true;
    *out = owner.release();
    return IRRED_OK;
  } catch (const std::bad_alloc &) {
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t
irred_expansion_result_view(const irred_expansion_result *r,
                            const irred_expansion_row **rows, uint64_t *count,
                            uint32_t *status, uint32_t *num,
                            uint64_t *callbacks, uint64_t *segments) {
  if (aligned(rows))
    *rows = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(status))
    *status = 0;
  if (aligned(num))
    *num = 0;
  if (aligned(callbacks))
    *callbacks = 0;
  if (aligned(segments))
    *segments = 0;
  if (!aligned(r) || !aligned(rows) || !aligned(count) || !aligned(status) ||
      !aligned(num) || !aligned(callbacks) || !aligned(segments))
    return IRRED_INVALID_INPUT;
  *rows = r->rows.data();
  *count = r->rows.size();
  *status = r->status;
  *num = r->numerical_status;
  *callbacks = r->callbacks;
  *segments = r->segments;
  return IRRED_OK;
}
extern "C" uint32_t
irred_expansion_result_models(const irred_expansion_result *r,
                              const irred_expansion_model_view **v,
                              uint64_t *count, uint32_t *available) {
  if (aligned(v))
    *v = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(available))
    *available = 0;
  if (!aligned(r) || !aligned(v) || !aligned(count) || !aligned(available))
    return IRRED_INVALID_INPUT;
  *v = r->models.data();
  *count = r->models.size();
  *available = r->source_available;
  return IRRED_OK;
}
extern "C" uint32_t
irred_expansion_result_queries(const irred_expansion_result *r,
                               const irred_expansion_request **v,
                               uint64_t *count, uint32_t *available) {
  if (aligned(v))
    *v = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(available))
    *available = 0;
  if (!aligned(r) || !aligned(v) || !aligned(count) || !aligned(available))
    return IRRED_INVALID_INPUT;
  *v = r->queries.data();
  *count = r->queries.size();
  *available = r->source_available;
  return IRRED_OK;
}
extern "C" uint32_t
irred_expansion_result_nodes(const irred_expansion_result *r,
                             const irred_expansion_node **v, uint64_t *count) {
  if (aligned(v))
    *v = nullptr;
  if (aligned(count))
    *count = 0;
  if (!aligned(r) || !aligned(v) || !aligned(count))
    return IRRED_INVALID_INPUT;
  *v = r->nodes.data();
  *count = r->nodes.size();
  return IRRED_OK;
}
extern "C" uint32_t irred_expansion_result_destroy(irred_expansion_result *r) {
  if (!r)
    return IRRED_OK;
  if (!aligned(r))
    return IRRED_INVALID_INPUT;
  delete r;
  return IRRED_OK;
}
