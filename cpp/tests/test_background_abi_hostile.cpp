// Bridge-independent comparison against the native P01 model, plus hostile
// structural/resource/allocation contracts. Independent mathematical P01
// fixtures live in test_background_hostile; this suite qualifies transport.
#include "irred/abi.h"
#include "irred/background.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
namespace {
long long fail_after = -1, live = 0;
int checks = 0;
} // namespace
void *operator new(std::size_t n) {
  if (fail_after == 0)
    throw std::bad_alloc();
  if (fail_after > 0)
    --fail_after;
  auto p = std::malloc(n ? n : 1);
  if (!p)
    throw std::bad_alloc();
  ++live;
  return p;
}
void operator delete(void *p) noexcept {
  if (p) {
    --live;
    std::free(p);
  }
}
void operator delete(void *p, std::size_t) noexcept { ::operator delete(p); }
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete[](void *p, std::size_t) noexcept { ::operator delete(p); }
namespace {
void check(bool ok, const char *why) {
  ++checks;
  if (!ok)
    throw std::runtime_error(why);
}
std::string_view text(cosmo_bytes b) {
  return {reinterpret_cast<const char *>(b.data),
          static_cast<std::size_t>(b.length)};
}
struct Fixture {
  std::array<cosmo_background_parameters, 2> parameters{
      {{COSMO_BACKGROUND_MODEL_FLAT_LCDM_LATE_V1, 0, 70, .3, 0},
       {COSMO_BACKGROUND_MODEL_CONSTANT_Q_FLAT_V1, 0, 70, 0, -1}}};
  std::array<cosmo_background_query, 3> queries{
      {{0, 0, COSMO_BACKGROUND_CONVENTION_GEOMETRIC_SAME_REDSHIFT, 0},
       {1, 1, COSMO_BACKGROUND_CONVENTION_GEOMETRIC_SAME_REDSHIFT, 0},
       {1, .9, COSMO_BACKGROUND_CONVENTION_RELEASED_ZHD_ZHEL, 0}}};
  cosmo_background_policy policy{sizeof(cosmo_background_policy),
                                 COSMO_ABI_VERSION,
                                 30,
                                 0,
                                 2,
                                 3,
                                 6,
                                 6 * sizeof(cosmo_background_slot),
                                 2000000,
                                 100000,
                                 1e-13,
                                 1e-12};
  cosmo_background_batch batch{sizeof(cosmo_background_batch),
                               COSMO_ABI_VERSION,
                               parameters.data(),
                               parameters.size(),
                               sizeof(parameters),
                               queries.data(),
                               queries.size(),
                               sizeof(queries)};
};
} // namespace
int main() {
  using namespace irred::cosmology;
  Fixture f;
  cosmo_background_result *result = nullptr;
  check(cosmo_background_evaluate(&f.batch, &f.policy, &result) == COSMO_OK &&
            result,
        "ordinary owned batch");
  const cosmo_background_slot *slots = nullptr;
  uint64_t count = 0;
  check(cosmo_background_result_view(result, &slots, &count) == COSMO_OK &&
            count == 6,
        "parameter major view");
  std::array<Query, 3> queries{{{0, 0, Convention::geometric_same_redshift},
                                {1, 1, Convention::geometric_same_redshift},
                                {1, .9, Convention::released_zhd_zhel}}};
  std::size_t remaining = f.policy.maximum_total_evaluations;
  for (std::size_t i = 0; i < 2; ++i) {
    const auto &p = f.parameters[i];
    auto model = prepare(
        {static_cast<Model>(p.model), p.h0_km_s_mpc, p.omega_m, p.constant_q});
    auto native = model.evaluate_batch(
        queries, {{1e-13, 1e-12, 100000, 30}, 3, remaining});
    for (std::size_t j = 0; j < 3; ++j) {
      const auto &a = slots[i * 3 + j];
      const auto &n = native.slots[j];
      remaining -= n.evaluations;
      check(a.parameter_index == i && a.query_index == j,
            "ordered pair indices");
      check(a.status == static_cast<uint32_t>(n.status) &&
                a.numerical_status == static_cast<uint32_t>(n.numerical_status),
            "direct status parity");
      check(a.evaluations == n.evaluations, "callback count parity");
      check(a.expansion_e == n.expansion_E && a.h_km_s_mpc == n.h_km_s_mpc &&
                a.radial_integral == n.radial_integral,
            "direct E H integral parity");
      check(a.radial_mpc == n.radial_mpc &&
                a.luminosity_mpc == n.luminosity_mpc &&
                a.angular_diameter_mpc == n.angular_diameter_mpc,
            "direct distances parity");
      check(a.dimensionless_luminosity_shape ==
                    n.dimensionless_luminosity_shape &&
                a.lookback_seconds == n.lookback_seconds &&
                a.volume_mpc3_per_sr_per_redshift ==
                    n.volume_mpc3_per_sr_per_redshift,
            "units/output parity");
      check(a.deceleration_q == n.deceleration_q && a.jerk == n.jerk,
            "kinematic parity");
      check(text(a.model_id) == model.model_id() &&
                text(a.constants_id) == Background::constants_id &&
                text(a.luminosity_equation_id) == n.luminosity_equation_id &&
                text(a.shape_equation_id) == n.shape_equation_id,
            "scientific identity parity");
    }
  }
  check(slots[0].radial_mpc == 0 && slots[0].luminosity_mpc == 0 &&
            slots[0].jerk == 1,
        "zero geometry keeps valid kinematics");
  check(text(slots[1].luminosity_equation_id) !=
            text(slots[2].luminosity_equation_id),
        "released luminosity has distinct equation ID");
  f.parameters[0].h0_km_s_mpc = 99;
  f.queries[0].z_expansion = 4;
  check(slots[0].parameters.h0_km_s_mpc == 70 &&
            slots[0].query.z_expansion == 0,
        "result owns source copies");
  check(cosmo_background_result_destroy(result) == COSMO_OK, "owned destroy");
  result = nullptr;
  f = Fixture{};
  // Fixture assignment requires resetting descriptor pointers to this object.
  f.batch.parameters = f.parameters.data();
  f.batch.queries = f.queries.data();
  auto reject = [&](cosmo_background_batch batch,
                    cosmo_background_policy policy, uint32_t expected) {
    result = reinterpret_cast<cosmo_background_result *>(1);
    auto code = cosmo_background_evaluate(&batch, &policy, &result);
    check(code == expected && !result, "structural rejection clears output");
  };
  auto bad = f.batch;
  bad.abi_version++;
  reject(bad, f.policy, COSMO_ABI_MISMATCH);
  auto policy = f.policy;
  policy.abi_version++;
  reject(f.batch, policy, COSMO_ABI_MISMATCH);
  bad = f.batch;
  bad.struct_size--;
  reject(bad, f.policy, COSMO_INVALID_INPUT);
  bad = f.batch;
  bad.parameter_byte_length--;
  reject(bad, f.policy, COSMO_INVALID_INPUT);
  bad = f.batch;
  bad.query_byte_length--;
  reject(bad, f.policy, COSMO_INVALID_INPUT);
  bad = f.batch;
  bad.parameters = nullptr;
  reject(bad, f.policy, COSMO_INVALID_INPUT);
  bad = f.batch;
  bad.queries = reinterpret_cast<const cosmo_background_query *>(1);
  reject(bad, f.policy, COSMO_INVALID_INPUT);
  bad = f.batch;
  bad.parameter_count = std::numeric_limits<uint64_t>::max();
  reject(bad, f.policy, COSMO_INVALID_INPUT);
  policy = f.policy;
  policy.maximum_slots = 5;
  reject(f.batch, policy, COSMO_INVALID_INPUT);
  policy = f.policy;
  policy.maximum_native_output_bytes = 6 * sizeof(cosmo_background_slot) - 1;
  reject(f.batch, policy, COSMO_INVALID_INPUT);
  policy = f.policy;
  policy.relative_tolerance = std::numeric_limits<double>::quiet_NaN();
  reject(f.batch, policy, COSMO_INVALID_INPUT);
  policy = f.policy;
  policy.maximum_depth = 0;
  reject(f.batch, policy, COSMO_INVALID_INPUT);
  f.parameters[0].reserved = 1;
  reject(f.batch, f.policy, COSMO_INVALID_INPUT);
  f.parameters[0].reserved = 0;
  f.queries[0].reserved = 1;
  reject(f.batch, f.policy, COSMO_INVALID_INPUT);
  f.queries[0].reserved = 0;
  policy = f.policy;
  policy.maximum_slots = 5;
  fail_after = 0;
  auto code = cosmo_background_evaluate(&f.batch, &policy, &result);
  fail_after = -1;
  check(code == COSMO_INVALID_INPUT && !result,
        "product cap precedes allocation");
  f.parameters[0].model = 99;
  check(cosmo_background_evaluate(&f.batch, &f.policy, &result) == COSMO_OK &&
            result,
        "unknown model is tagged native failure");
  check(cosmo_background_result_view(result, &slots, &count) == COSMO_OK &&
            slots[0].status == COSMO_BACKGROUND_STATUS_INVALID_INPUT &&
            slots[0].expansion_e == 0,
        "failed slot has no physical payload");
  cosmo_background_result_destroy(result);
  result = nullptr;
  f.parameters[0].model = COSMO_BACKGROUND_MODEL_FLAT_LCDM_LATE_V1;
  policy = f.policy;
  policy.maximum_total_evaluations = 1;
  check(cosmo_background_evaluate(&f.batch, &policy, &result) == COSMO_OK,
        "work failure stays in owned slots");
  cosmo_background_result_view(result, &slots, &count);
  uint64_t used = 0;
  for (std::size_t i = 0; i < count; ++i) {
    used += slots[i].evaluations;
    if (i % 3)
      check(slots[i].status == COSMO_BACKGROUND_STATUS_WORK_LIMIT &&
                slots[i].numerical_status ==
                    COSMO_NUMERICAL_STATUS_WORK_LIMIT &&
                slots[i].radial_mpc == 0,
            "global work failure cause/payload");
  }
  check(used <= 1 && slots[0].status == COSMO_BACKGROUND_STATUS_OK &&
            slots[3].status == COSMO_BACKGROUND_STATUS_OK,
        "zero queries require no callback budget");
  cosmo_background_result_destroy(result);
  result = nullptr;
  bool succeeded = false;
  int failures = 0;
  for (int point = 0; point < 100; ++point) {
    const auto before = live;
    fail_after = point;
    code = cosmo_background_evaluate(&f.batch, &f.policy, &result);
    fail_after = -1;
    if (code == COSMO_OK) {
      check(result, "successful allocation sweep result");
      cosmo_background_result_destroy(result);
      result = nullptr;
      check(live == before, "success cleanup balanced");
      succeeded = true;
      break;
    }
    check(code == COSMO_ALLOCATION_FAILURE && !result,
          "bad_alloc cannot cross ABI");
    check(live == before, "partial allocation cleanup balanced");
    ++failures;
  }
  check(succeeded && failures >= 5,
        "allocation sweep reaches every path then success");
  check(cosmo_background_result_destroy(nullptr) == COSMO_OK, "null destroy");
  check(cosmo_background_result_view(nullptr, &slots, &count) ==
                COSMO_INVALID_INPUT &&
            !slots && !count,
        "invalid view initializes outputs");
  std::printf("{\"suite\":\"background_abi_hostile\",\"checks\":%d,"
              "\"allocation_failure_points\":%d,\"passed\":true}\n",
              checks, failures);
}
