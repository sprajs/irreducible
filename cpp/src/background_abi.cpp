#include "irred/abi.h"
#include "irred/background.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <type_traits>
#include <vector>
namespace {
namespace physics = irred::cosmology;
static_assert(static_cast<uint32_t>(physics::Model::flat_cpl_late_v1) ==
              COSMO_BACKGROUND_V2_MODEL_FLAT_CPL_LATE_V1);
static_assert(sizeof(cosmo_background_parameters) == 32);
static_assert(sizeof(cosmo_background_parameters_v2) == 48);
static_assert(static_cast<uint32_t>(physics::Model::flat_lcdm_late_v1) ==
              COSMO_BACKGROUND_MODEL_FLAT_LCDM_LATE_V1);
static_assert(static_cast<uint32_t>(physics::Model::constant_q_flat_v1) ==
              COSMO_BACKGROUND_MODEL_CONSTANT_Q_FLAT_V1);
static_assert(
    static_cast<uint32_t>(physics::Convention::geometric_same_redshift) ==
    COSMO_BACKGROUND_CONVENTION_GEOMETRIC_SAME_REDSHIFT);
static_assert(static_cast<uint32_t>(physics::Convention::released_zhd_zhel) ==
              COSMO_BACKGROUND_CONVENTION_RELEASED_ZHD_ZHEL);
static_assert(static_cast<uint32_t>(physics::Status::ok) ==
              COSMO_BACKGROUND_STATUS_OK);
static_assert(static_cast<uint32_t>(physics::Status::invalid_input) ==
              COSMO_BACKGROUND_STATUS_INVALID_INPUT);
static_assert(static_cast<uint32_t>(physics::Status::unsupported_domain) ==
              COSMO_BACKGROUND_STATUS_UNSUPPORTED_DOMAIN);
static_assert(static_cast<uint32_t>(physics::Status::incompatible_convention) ==
              COSMO_BACKGROUND_STATUS_INCOMPATIBLE_CONVENTION);
static_assert(static_cast<uint32_t>(physics::Status::numerical_failure) ==
              COSMO_BACKGROUND_STATUS_NUMERICAL_FAILURE);
static_assert(static_cast<uint32_t>(physics::Status::work_limit) ==
              COSMO_BACKGROUND_STATUS_WORK_LIMIT);
template <class T> bool aligned(const T *p) {
  return p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0;
}
template <class T>
bool bounded(const T *p, uint64_t n, uint64_t bytes, uint64_t cap) {
  return n <= cap && n <= std::numeric_limits<size_t>::max() / sizeof(T) &&
         bytes == n * sizeof(T) &&
         (!n || (p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0));
}
cosmo_bytes view(std::string_view text) {
  return {reinterpret_cast<const uint8_t *>(text.data()), text.size()};
}
bool valid(const cosmo_background_policy &p) {
  return p.struct_size == sizeof(p) && p.abi_version == COSMO_ABI_VERSION &&
         !p.reserved && p.maximum_parameters <= 4096 &&
         p.maximum_queries <= 4096 &&
         p.maximum_slots <= COSMO_MAX_BATCH_ELEMENTS &&
         p.maximum_native_output_bytes <= 512u * 1024u * 1024u &&
         p.maximum_total_evaluations <= std::numeric_limits<size_t>::max() &&
         p.maximum_evaluations_per_integral >= 3 &&
         p.maximum_evaluations_per_integral <=
             std::numeric_limits<size_t>::max() &&
         p.maximum_depth > 0 && p.maximum_depth <= 64 &&
         std::isfinite(p.absolute_tolerance) &&
         std::isfinite(p.relative_tolerance) && p.absolute_tolerance >= 0 &&
         p.relative_tolerance >= 0 &&
         (p.absolute_tolerance > 0 || p.relative_tolerance > 0);
}
} // namespace
struct cosmo_background_result_v2 {
  std::vector<cosmo_background_slot_v2> slots;
};
struct cosmo_background_result {
  std::vector<cosmo_background_slot> slots;
};
template <class Batch, class Result, class WireSlot>
uint32_t evaluate_background(const Batch *b, const cosmo_background_policy *p,
                             Result **out) {
  if (aligned(out))
    *out = nullptr;
  if (!aligned(out) || !aligned(b) || !aligned(p))
    return COSMO_INVALID_INPUT;
  if (b->abi_version != COSMO_ABI_VERSION ||
      p->abi_version != COSMO_ABI_VERSION)
    return COSMO_ABI_MISMATCH;
  if (b->struct_size != sizeof(*b) || !valid(*p) ||
      !bounded(b->parameters, b->parameter_count, b->parameter_byte_length,
               p->maximum_parameters) ||
      !bounded(b->queries, b->query_count, b->query_byte_length,
               p->maximum_queries))
    return COSMO_INVALID_INPUT;
  // Product and output byte limits precede input copies, result allocation and
  // physics.
  if (b->query_count && b->parameter_count > p->maximum_slots / b->query_count)
    return COSMO_INVALID_INPUT;
  const auto total = b->parameter_count * b->query_count;
  if (total > p->maximum_native_output_bytes / sizeof(WireSlot) ||
      total > std::numeric_limits<size_t>::max() / sizeof(WireSlot))
    return COSMO_INVALID_INPUT;
  // Also guard the sum of owned output, copied queries and one temporary native
  // per-model batch. Count caps bound these work buffers independently.
  const auto output_bytes = static_cast<size_t>(total) * sizeof(WireSlot);
  if (b->query_count > (std::numeric_limits<size_t>::max() - output_bytes) /
                           sizeof(physics::Query))
    return COSMO_INVALID_INPUT;
  const auto with_queries = output_bytes + static_cast<size_t>(b->query_count) *
                                               sizeof(physics::Query);
  if (b->query_count > (std::numeric_limits<size_t>::max() - with_queries) /
                           sizeof(physics::Slot))
    return COSMO_INVALID_INPUT;
  for (uint64_t i = 0; i < b->parameter_count; ++i)
    if (b->parameters[i].reserved)
      return COSMO_INVALID_INPUT;
  for (uint64_t i = 0; i < b->query_count; ++i)
    if (b->queries[i].reserved)
      return COSMO_INVALID_INPUT;
  try {
    std::vector<physics::Query> queries;
    queries.reserve(b->query_count);
    for (uint64_t i = 0; i < b->query_count; ++i)
      queries.push_back(
          {b->queries[i].z_expansion, b->queries[i].z_observer,
           static_cast<physics::Convention>(b->queries[i].convention)});
    auto result = std::make_unique<Result>();
    result->slots.resize(total);
    size_t remaining = p->maximum_total_evaluations;
    for (uint64_t i = 0; i < b->parameter_count; ++i) {
      const auto &source = b->parameters[i];
      // v1 carries no CPL parameters; unknown IDs stay unsupported per row.
      auto model = [&] {
        if constexpr (std::is_same_v<Batch, cosmo_background_batch_v2>) {
          if (source.model == COSMO_BACKGROUND_V2_MODEL_FLAT_CPL_LATE_V1 &&
              source.constant_q == 0)
            return physics::prepare_cpl(
                {source.h0_km_s_mpc, source.omega_m, source.w0, source.wa});
          auto value = physics::Parameters{
              static_cast<physics::Model>(source.model), source.h0_km_s_mpc,
              source.omega_m, source.constant_q};
          value.w0 = source.w0;
          value.wa = source.wa;
          return physics::prepare(value);
        } else {
          const auto model_id = source.model <= 1 ? source.model : UINT32_MAX;
          return physics::prepare({static_cast<physics::Model>(model_id),
                                   source.h0_km_s_mpc, source.omega_m,
                                   source.constant_q});
        }
      }();
      physics::Policy policy{
          {p->absolute_tolerance, p->relative_tolerance,
           static_cast<size_t>(p->maximum_evaluations_per_integral),
           p->maximum_depth},
          static_cast<size_t>(p->maximum_queries),
          remaining};
      auto batch = model.evaluate_batch(queries, policy);
      if (batch.status == physics::Status::ok &&
          batch.slots.size() != queries.size())
        return COSMO_EXCEPTION;
      for (uint64_t j = 0; j < b->query_count; ++j) {
        auto &slot = result->slots[i * b->query_count + j];
        slot.parameter_index = i;
        slot.query_index = j;
        slot.parameters = source;
        slot.query = b->queries[j];
        slot.model_id = view(model.model_id());
        slot.constants_id = view(physics::Background::constants_id);
        slot.radial_equation_id = view(physics::Background::radial_equation_id);
        slot.numerical_status = COSMO_NUMERICAL_STATUS_INVALID_INPUT;
        if (batch.status != physics::Status::ok) {
          slot.status = static_cast<uint32_t>(batch.status);
          continue;
        }
        const auto &value = batch.slots[j];
        if (value.evaluations > remaining)
          return COSMO_EXCEPTION;
        remaining -= value.evaluations;
        slot.status = static_cast<uint32_t>(value.status);
        slot.numerical_status = static_cast<uint32_t>(value.numerical_status);
        slot.luminosity_equation_id = view(value.luminosity_equation_id);
        slot.shape_equation_id = view(value.shape_equation_id);
        slot.evaluations = value.evaluations;
        // Only successful scientific bundles expose finite physical payloads.
        if (value.status != physics::Status::ok)
          continue;
        slot.expansion_e = value.expansion_E;
        slot.h_km_s_mpc = value.h_km_s_mpc;
        slot.radial_integral = value.radial_integral;
        slot.radial_mpc = value.radial_mpc;
        slot.transverse_mpc = value.transverse_mpc;
        slot.angular_diameter_mpc = value.angular_diameter_mpc;
        slot.luminosity_mpc = value.luminosity_mpc;
        slot.dimensionless_luminosity_shape =
            value.dimensionless_luminosity_shape;
        slot.lookback_seconds = value.lookback_seconds;
        slot.volume_mpc3_per_sr_per_redshift =
            value.volume_mpc3_per_sr_per_redshift;
        slot.deceleration_q = value.deceleration_q;
        slot.jerk = value.jerk;
        slot.radial_integral_error = value.radial_integral_error;
        slot.lookback_integral_error = value.lookback_integral_error;
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
extern "C" uint32_t cosmo_background_evaluate(const cosmo_background_batch *b,
                                              const cosmo_background_policy *p,
                                              cosmo_background_result **out) {
  return evaluate_background<cosmo_background_batch, cosmo_background_result,
                             cosmo_background_slot>(b, p, out);
}
extern "C" uint32_t
cosmo_background_evaluate_v2(const cosmo_background_batch_v2 *b,
                             const cosmo_background_policy *p,
                             cosmo_background_result_v2 **out) {
  return evaluate_background<cosmo_background_batch_v2,
                             cosmo_background_result_v2,
                             cosmo_background_slot_v2>(b, p, out);
}
extern "C" uint32_t
cosmo_background_result_view(const cosmo_background_result *r,
                             const cosmo_background_slot **slots,
                             uint64_t *length) {
  if (aligned(slots))
    *slots = nullptr;
  if (aligned(length))
    *length = 0;
  if (!aligned(r) || !aligned(slots) || !aligned(length))
    return COSMO_INVALID_INPUT;
  *slots = r->slots.data();
  *length = r->slots.size();
  return COSMO_OK;
}
extern "C" uint32_t
cosmo_background_result_destroy(cosmo_background_result *r) {
  delete r;
  return COSMO_OK;
}

extern "C" uint32_t
cosmo_background_result_v2_view(const cosmo_background_result_v2 *r,
                                const cosmo_background_slot_v2 **slots,
                                uint64_t *length) {
  if (aligned(slots))
    *slots = nullptr;
  if (aligned(length))
    *length = 0;
  if (!aligned(r) || !aligned(slots) || !aligned(length))
    return COSMO_INVALID_INPUT;
  *slots = r->slots.data();
  *length = r->slots.size();
  return COSMO_OK;
}
extern "C" uint32_t
cosmo_background_result_v2_destroy(cosmo_background_result_v2 *r) {
  delete r;
  return COSMO_OK;
}
