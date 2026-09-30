#include "irred/abi.h"
#include "irred/piecewise_background.hpp"
#include <algorithm>
#include <limits>
#include <memory>
#include <new>
#include <string_view>
#include <vector>

namespace {
namespace physics = irred::cosmology;
static_assert(sizeof(cosmo_piecewise_parameters) == 48);
static_assert(sizeof(cosmo_piecewise_batch) == 56);
static_assert(sizeof(cosmo_piecewise_policy) == 56);
static_assert(sizeof(cosmo_piecewise_slot) == 328);
static_assert(static_cast<uint32_t>(physics::QConvention::not_assessed) ==
              COSMO_PIECEWISE_Q_CONVENTION_NOT_ASSESSED);
static_assert(
    static_cast<uint32_t>(physics::QConvention::interior_constant_bin) ==
    COSMO_PIECEWISE_Q_CONVENTION_INTERIOR_CONSTANT_BIN);
static_assert(
    static_cast<uint32_t>(physics::QConvention::right_limit_at_internal_jump) ==
    COSMO_PIECEWISE_Q_CONVENTION_RIGHT_LIMIT_AT_INTERNAL_JUMP);
static_assert(
    static_cast<uint32_t>(physics::QConvention::right_limit_at_zero) ==
    COSMO_PIECEWISE_Q_CONVENTION_RIGHT_LIMIT_AT_ZERO);
static_assert(
    static_cast<uint32_t>(physics::QConvention::left_limit_at_final_endpoint) ==
    COSMO_PIECEWISE_Q_CONVENTION_LEFT_LIMIT_AT_FINAL_ENDPOINT);
static_assert(static_cast<uint32_t>(physics::JerkAvailability::not_assessed) ==
              COSMO_PIECEWISE_JERK_AVAILABILITY_NOT_ASSESSED);
static_assert(
    static_cast<uint32_t>(physics::JerkAvailability::ordinary_within_bin) ==
    COSMO_PIECEWISE_JERK_AVAILABILITY_ORDINARY_WITHIN_BIN);
static_assert(
    static_cast<uint32_t>(physics::JerkAvailability::one_sided_endpoint) ==
    COSMO_PIECEWISE_JERK_AVAILABILITY_ONE_SIDED_ENDPOINT);
static_assert(
    static_cast<uint32_t>(physics::JerkAvailability::unavailable_at_jump) ==
    COSMO_PIECEWISE_JERK_AVAILABILITY_UNAVAILABLE_AT_JUMP);
template <class T> bool aligned(const T *p) {
  return p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0;
}
template <class T> bool descriptor(const T *p, uint64_t count, uint64_t bytes) {
  return count <= std::numeric_limits<size_t>::max() / sizeof(T) &&
         bytes == count * sizeof(T) && (!count || aligned(p));
}
cosmo_bytes view(std::string_view s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
bool add_bytes(size_t &total, uint64_t count, size_t width) {
  if (count > (std::numeric_limits<size_t>::max() - total) / width)
    return false;
  total += static_cast<size_t>(count) * width;
  return true;
}
} // namespace
struct cosmo_piecewise_result {
  std::vector<cosmo_piecewise_slot> slots;
  uint32_t status = COSMO_BACKGROUND_STATUS_OK;
  uint32_t numerical_status = COSMO_NUMERICAL_STATUS_OK;
  uint64_t segments_processed = 0;
};
extern "C" uint32_t cosmo_piecewise_evaluate(const cosmo_piecewise_batch *b,
                                             const cosmo_piecewise_policy *p,
                                             cosmo_piecewise_result **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(b) || !aligned(p))
    return COSMO_INVALID_INPUT;
  if (b->abi_version != COSMO_ABI_VERSION ||
      p->abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (b->struct_size != sizeof(*b) || p->struct_size != sizeof(*p) ||
      p->reserved || p->maximum_parameters > 4096 ||
      p->maximum_queries > 4096 || p->maximum_slots > 1048576 ||
      p->maximum_native_output_bytes > 512u * 1024u * 1024u ||
      p->maximum_total_segment_visits > std::numeric_limits<size_t>::max() ||
      !descriptor(b->parameters, b->parameter_count, b->parameter_bytes) ||
      !descriptor(b->queries, b->query_count, b->query_bytes))
    return COSMO_INVALID_INPUT;
  // Cap checks precede dereferencing/copying any row. Valid resource exhaustion
  // has an independently owned empty failure; malformed storage has null out.
  bool capped = b->parameter_count > p->maximum_parameters ||
                b->query_count > p->maximum_queries;
  uint64_t total = 0;
  size_t bytes = sizeof(cosmo_piecewise_result);
  if (!capped) {
    if (b->query_count &&
        b->parameter_count > p->maximum_slots / b->query_count)
      capped = true;
    else {
      total = b->parameter_count * b->query_count;
      if (!add_bytes(bytes, total, sizeof(cosmo_piecewise_slot)) ||
          !add_bytes(bytes, b->parameter_count,
                     sizeof(cosmo_piecewise_parameters)) ||
          !add_bytes(bytes, b->query_count, sizeof(physics::Query)) ||
          !add_bytes(bytes, b->query_count, sizeof(physics::PiecewiseSlot)))
        return COSMO_INVALID_INPUT;
      capped = bytes > p->maximum_native_output_bytes;
    }
  }
  if (!capped)
    for (uint64_t j = 0; j < b->query_count; ++j)
      if (b->queries[j].reserved)
        return COSMO_INVALID_INPUT;
  try {
    auto result = std::make_unique<cosmo_piecewise_result>();
    if (capped) {
      result->status = COSMO_BACKGROUND_STATUS_WORK_LIMIT;
      result->numerical_status = COSMO_NUMERICAL_STATUS_WORK_LIMIT;
      *out = result.release();
      return COSMO_OK;
    }
    std::vector<cosmo_piecewise_parameters> parameters;
    if (b->parameter_count)
      parameters.assign(b->parameters, b->parameters + b->parameter_count);
    std::vector<physics::Query> queries;
    queries.reserve(b->query_count);
    for (uint64_t j = 0; j < b->query_count; ++j)
      queries.push_back(
          {b->queries[j].z_expansion, b->queries[j].z_observer,
           static_cast<physics::Convention>(b->queries[j].convention)});
    result->slots.resize(total);
    size_t remaining = static_cast<size_t>(p->maximum_total_segment_visits);
    for (uint64_t i = 0; i < b->parameter_count; ++i) {
      const auto &source = parameters[i];
      std::array<double, 5> q;
      std::copy_n(source.q, 5, q.begin());
      auto model = physics::prepare_piecewise_q({source.h0_km_s_mpc, q});
      auto batch = model.evaluate_batch(
          queries, {static_cast<size_t>(p->maximum_queries), remaining});
      if (batch.segments_processed > remaining)
        return COSMO_EXCEPTION;
      remaining -= batch.segments_processed;
      result->segments_processed += batch.segments_processed;
      if (batch.status == physics::Status::ok &&
          batch.slots.size() != queries.size())
        return COSMO_EXCEPTION;
      for (uint64_t j = 0; j < b->query_count; ++j) {
        auto &row = result->slots[i * b->query_count + j];
        row.struct_size = sizeof(row);
        row.abi_version = COSMO_ABI_VERSION;
        row.parameter_index = i;
        row.query_index = j;
        row.parameters = source;
        row.query = b->queries[j];
        row.status = static_cast<uint32_t>(batch.status);
        row.numerical_status = COSMO_NUMERICAL_STATUS_INVALID_INPUT;
        row.model_id = view(physics::PiecewiseBackground::model_id);
        row.constant_set_id = view(physics::PiecewiseBackground::constants_id);
        row.radial_equation_id = view(physics::Background::radial_equation_id);
        if (row.query.convention ==
            COSMO_BACKGROUND_CONVENTION_GEOMETRIC_SAME_REDSHIFT) {
          row.luminosity_equation_id =
              view("P01/flat-geometric-luminosity-distance/v1");
          row.shape_equation_id =
              view("P01/geometric-same-redshift-dimensionless-shape/v1");
        } else if (row.query.convention ==
                   COSMO_BACKGROUND_CONVENTION_RELEASED_ZHD_ZHEL) {
          row.luminosity_equation_id =
              view("P01/released-zhd-zhel-luminosity-prefactor/v1");
          row.shape_equation_id =
              view("P01/released-zhd-zhel-intercept-free-shape/v1");
        }
        if (batch.status != physics::Status::ok) {
          if (batch.status == physics::Status::work_limit)
            row.numerical_status = COSMO_NUMERICAL_STATUS_WORK_LIMIT;
          continue;
        }
        const auto &v = batch.slots[j];
        row.status = static_cast<uint32_t>(v.status);
        row.numerical_status = static_cast<uint32_t>(v.numerical_status);
        row.segments_processed = v.segments_processed;
        if (v.status != physics::Status::ok)
          continue;
        row.has_geometry = 1;
        row.bin = v.bin;
        row.q_convention = static_cast<uint32_t>(v.q_convention);
        row.jerk_availability = static_cast<uint32_t>(v.jerk_availability);
        if (v.assigned_q) {
          row.has_q = 1;
          row.assigned_q = *v.assigned_q;
        }
        if (v.jerk) {
          row.has_jerk = 1;
          row.jerk = *v.jerk;
        }
        if (v.q0_within_piecewise_model) {
          row.has_q0 = 1;
          row.q0_within_piecewise_model = *v.q0_within_piecewise_model;
        }
        const auto &g = v.geometry;
        row.expansion_E = g.expansion_E;
        row.h_km_s_mpc = g.h_km_s_mpc;
        row.radial_integral = g.radial_integral;
        row.radial_mpc = g.radial_mpc;
        row.transverse_mpc = g.transverse_mpc;
        row.angular_diameter_mpc = g.angular_diameter_mpc;
        row.luminosity_mpc = g.luminosity_mpc;
        row.dimensionless_luminosity_shape = g.dimensionless_luminosity_shape;
        row.lookback_seconds = g.lookback_seconds;
        row.volume_mpc3_per_sr_per_redshift = g.volume_mpc3_per_sr_per_redshift;
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
extern "C" uint32_t
cosmo_piecewise_result_view(const cosmo_piecewise_result *r,
                            const cosmo_piecewise_slot **rows, uint64_t *count,
                            uint32_t *status, uint32_t *numerical_status,
                            uint64_t *segments) {
  if (aligned(rows))
    *rows = nullptr;
  if (aligned(count))
    *count = 0;
  if (aligned(status))
    *status = COSMO_BACKGROUND_STATUS_INVALID_INPUT;
  if (aligned(numerical_status))
    *numerical_status = COSMO_NUMERICAL_STATUS_INVALID_INPUT;
  if (aligned(segments))
    *segments = 0;
  if (!aligned(r) || !aligned(rows) || !aligned(count) || !aligned(status) ||
      !aligned(numerical_status) || !aligned(segments))
    return COSMO_INVALID_INPUT;
  *rows = r->slots.data();
  *count = r->slots.size();
  *status = r->status;
  *numerical_status = r->numerical_status;
  *segments = r->segments_processed;
  return COSMO_OK;
}
extern "C" uint32_t cosmo_piecewise_result_destroy(cosmo_piecewise_result *r) {
  delete r;
  return COSMO_OK;
}
