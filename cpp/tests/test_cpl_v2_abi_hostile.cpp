// Independent v2 transport review. Fixture provenance remains synthetic;
// native physics independent comparisons reside in CPL peer suites.
// Transport adversaries by consumer author, independently of bridge
// implementation. Generated release-shaped bytes test interface semantics, not
// actual survey science.
#include "irred/abi.h"
#include "irred/background.hpp"
#include "irred/observations.hpp"
#include "irred/supernova.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
static std::int64_t fail_after = -1, live = 0;
static int checks = 0;
void *operator new(std::size_t n) {
  if (fail_after == 0)
    throw std::bad_alloc();
  if (fail_after > 0)
    --fail_after;
  void *p = std::malloc(n ? n : 1);
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
#define CHECK(...)                                                             \
  do {                                                                         \
    ++checks;                                                                  \
    if (!(__VA_ARGS__))                                                        \
      throw std::runtime_error("Supernova ABI: " #__VA_ARGS__);                \
  } while (false)
using namespace irred::observations;
cosmo_bytes bytes(const std::string &s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
cosmo_f64_buffer doubles(const double *p, std::size_t n) {
  return {sizeof(cosmo_f64_buffer), COSMO_ABI_VERSION, 2, 0, p, n,
          n * sizeof(double)};
}
struct ObservationFixture {
  std::string hash = std::string(64, 'a'), covhash = std::string(64, 'b'),
              order = "synthetic explicit axes", a = "a", b = "b",
              event = "same", latent = "latent-one";
  cosmo_bytes ids[2] = {bytes(a), bytes(b)},
              events[2] = {bytes(event), bytes(event)};
  double values[2] = {2, -3}, matrix[4] = {4, 1, 1, 9}, z[2] = {.1, .2};
  uint8_t missing[2] = {0, 0};
  uint64_t quality[2] = {0, 0};
  cosmo_observation_policy op{2, 4, 4096};
  cosmo_observation_descriptor d{};
  ObservationFixture() {
    d.struct_size = sizeof(d);
    d.abi_version = COSMO_ABI_VERSION;
    d.profile = static_cast<uint32_t>(Profile::pantheon_plus_released_v1);
    d.role = static_cast<uint32_t>(Role::released_fitted_summary);
    d.unit = static_cast<uint32_t>(Unit::magnitude);
    d.calibration = static_cast<uint32_t>(Calibration::unknown);
    d.uncertainty = static_cast<uint32_t>(Uncertainty::covariance);
    d.uncertainty_unit =
        static_cast<uint32_t>(UncertaintyUnit::magnitude_squared);
    d.table_sha256 = bytes(hash);
    d.uncertainty_sha256 = bytes(covhash);
    d.ordering_provenance = bytes(order);
    d.measurement_ids = {ids, 2, sizeof ids};
    d.event_ids = {events, 2, sizeof events};
    d.uncertainty_axis_ids = d.measurement_ids;
    d.values = doubles(values, 2);
    d.zhd = d.zcmb = d.zhel = doubles(z, 2);
    d.zhd_missing = d.zcmb_missing = d.zhel_missing = {missing, 2, 2};
    d.uncertainty_matrix = doubles(matrix, 4);
    d.missing = {missing, 2, 2};
    d.quality = {quality, 2, sizeof quality};
  }
};
void background_test() {
  static_assert(sizeof(cosmo_background_parameters_v2) == 48);
  static_assert(offsetof(cosmo_background_parameters_v2, w0) == 32);
  cosmo_background_parameters_v2 parameters[] = {{2, 0, 70, .3, 0, -.9, .4},
                                                 {0, 0, 70, .3, 0, -1, 0},
                                                 {2, 0, 70, .3, 1, -.9, .4}};
  cosmo_background_query queries[] = {
      {0, 0, 0, 0}, {1, 1, 0, 0}, {2, 1.9, 1, 0}};
  cosmo_background_policy policy{sizeof(policy),
                                 COSMO_ABI_VERSION,
                                 30,
                                 0,
                                 3,
                                 3,
                                 9,
                                 9 * sizeof(cosmo_background_slot_v2),
                                 2000000,
                                 100000,
                                 1e-13,
                                 1e-12};
  cosmo_background_batch_v2 batch{
      sizeof(batch),      COSMO_ABI_VERSION, parameters, 3,
      sizeof(parameters), queries,           3,          sizeof(queries)};
  cosmo_background_result_v2 *result = nullptr;
  CHECK(cosmo_background_evaluate_v2(&batch, &policy, &result) == COSMO_OK &&
        result);
  const cosmo_background_slot_v2 *slots = nullptr;
  uint64_t count = 0;
  CHECK(cosmo_background_result_v2_view(result, &slots, &count) == COSMO_OK &&
        count == 9);
  std::array<irred::cosmology::Query, 3> native_queries{
      {{0, 0, irred::cosmology::Convention::geometric_same_redshift},
       {1, 1, irred::cosmology::Convention::geometric_same_redshift},
       {2, 1.9, irred::cosmology::Convention::released_zhd_zhel}}};
  std::size_t remaining = policy.maximum_total_evaluations;
  for (std::size_t m = 0; m < 2; ++m) {
    const auto &p = parameters[m];
    auto background = m == 0 ? irred::cosmology::prepare_cpl(
                                   {p.h0_km_s_mpc, p.omega_m, p.w0, p.wa})
                             : irred::cosmology::prepare(
                                   {irred::cosmology::Model::flat_lcdm_late_v1,
                                    p.h0_km_s_mpc, p.omega_m, p.constant_q});
    auto native = background.evaluate_batch(
        native_queries, {{1e-13, 1e-12, 100000, 30}, 3, remaining});
    for (std::size_t q = 0; q < 3; ++q) {
      const auto &a = slots[m * 3 + q];
      const auto &n = native.slots[q];
      remaining -= n.evaluations;
      CHECK(a.parameter_index == m && a.query_index == q &&
            a.parameters.w0 == p.w0 && a.parameters.wa == p.wa);
      CHECK(a.status == static_cast<uint32_t>(n.status) &&
            a.numerical_status == static_cast<uint32_t>(n.numerical_status));
      CHECK(a.radial_integral == n.radial_integral &&
            a.luminosity_mpc == n.luminosity_mpc &&
            a.lookback_seconds == n.lookback_seconds);
      CHECK(a.deceleration_q == n.deceleration_q && a.jerk == n.jerk &&
            a.evaluations == n.evaluations);
    }
  }
  CHECK(slots[6].status != 0 && slots[6].parameters.constant_q == 1 &&
        slots[6].parameters.w0 == -.9 && slots[6].radial_mpc == 0 &&
        slots[6].evaluations == 0);
  parameters[0].wa = -2;
  queries[0].z_expansion = 5;
  CHECK(slots[0].parameters.wa == .4 && slots[0].query.z_expansion == 0);
  CHECK(cosmo_background_result_v2_destroy(result) == COSMO_OK);
  parameters[0].wa = .4;
  queries[0].z_expansion = 0;
  auto before = live;
  int failures = 0;
  bool reached = false;
  for (int k = 0; k < 100; ++k) {
    result = nullptr;
    fail_after = k;
    auto status = cosmo_background_evaluate_v2(&batch, &policy, &result);
    fail_after = -1;
    if (status == COSMO_OK) {
      CHECK(result);
      cosmo_background_result_v2_destroy(result);
      CHECK(live == before);
      reached = true;
      break;
    }
    CHECK(status == COSMO_ALLOCATION_FAILURE && !result && live == before);
    ++failures;
  }
  CHECK(reached && failures > 0);
  std::printf("background_v2 allocation_failures %d\n", failures);
  auto impossible = batch;
  impossible.parameters = nullptr;
  impossible.parameter_count = UINT64_MAX;
  impossible.parameter_byte_length = 0;
  fail_after = 0;
  auto bad = cosmo_background_evaluate_v2(&impossible, &policy, &result);
  fail_after = -1;
  CHECK(bad == COSMO_INVALID_INPUT && !result && live == before);
  auto wrong = batch;
  wrong.abi_version++;
  CHECK(cosmo_background_evaluate_v2(&wrong, &policy, &result) ==
            COSMO_ABI_MISMATCH &&
        !result);
  wrong = batch;
  wrong.parameter_byte_length--;
  CHECK(cosmo_background_evaluate_v2(&wrong, &policy, &result) ==
            COSMO_INVALID_INPUT &&
        !result);
  wrong = batch;
  wrong.parameters = reinterpret_cast<const cosmo_background_parameters_v2 *>(
      reinterpret_cast<const char *>(parameters) + 1);
  CHECK(cosmo_background_evaluate_v2(&wrong, &policy, &result) ==
            COSMO_INVALID_INPUT &&
        !result);
  auto cap = policy;
  cap.maximum_native_output_bytes = 9 * sizeof(cosmo_background_slot_v2) - 1;
  CHECK(cosmo_background_evaluate_v2(&batch, &cap, &result) ==
            COSMO_INVALID_INPUT &&
        !result);
  cap = policy;
  cap.maximum_total_evaluations = 3;
  CHECK(cosmo_background_evaluate_v2(&batch, &cap, &result) == COSMO_OK &&
        result);
  CHECK(cosmo_background_result_v2_view(result, &slots, &count) == COSMO_OK);
  uint64_t calls = 0;
  for (uint64_t i = 0; i < count; ++i)
    calls += slots[i].evaluations;
  CHECK(calls <= 3 && slots[1].status != 0 && slots[4].status != 0);
  cosmo_background_result_v2_destroy(result);
  parameters[0].reserved = 1;
  CHECK(cosmo_background_evaluate_v2(&batch, &policy, &result) ==
            COSMO_INVALID_INPUT &&
        !result);
  parameters[0].reserved = 0;
  parameters[1].w0 = -.8;
  CHECK(cosmo_background_evaluate_v2(&batch, &policy, &result) == COSMO_OK &&
        result);
  CHECK(cosmo_background_result_v2_view(result, &slots, &count) == COSMO_OK &&
        slots[3].status != 0 && slots[3].parameters.w0 == -.8 &&
        slots[3].radial_mpc == 0);
  cosmo_background_result_v2_destroy(result);
  parameters[1].w0 = -1;
  parameters[0].w0 = std::numeric_limits<double>::quiet_NaN();
  CHECK(cosmo_background_evaluate_v2(&batch, &policy, &result) == COSMO_OK &&
        result);
  CHECK(cosmo_background_result_v2_view(result, &slots, &count) == COSMO_OK &&
        slots[0].status != 0 && std::isnan(slots[0].parameters.w0) &&
        slots[0].radial_mpc == 0);
  cosmo_background_result_v2_destroy(result);
  parameters[0].w0 = -.9;
  CHECK(cosmo_background_evaluate_v2(&batch, &policy, nullptr) ==
        COSMO_INVALID_INPUT);
  cosmo_background_parameters legacy{2, 0, 70, .3, 0};
  cosmo_background_batch v1{sizeof(v1), COSMO_ABI_VERSION, &legacy,
                            1,          sizeof(legacy),    queries,
                            1,          sizeof(queries[0])};
  cosmo_background_result *old = nullptr;
  CHECK(cosmo_background_evaluate(&v1, &policy, &old) == COSMO_OK && old);
  const cosmo_background_slot *os = nullptr;
  CHECK(cosmo_background_result_view(old, &os, &count) == COSMO_OK &&
        os[0].status != 0);
  cosmo_background_result_destroy(old);
}
void supernova_test() {
  static_assert(sizeof(cosmo_supernova_model_v2) == 40);
  static_assert(offsetof(cosmo_supernova_model_v2, w0) == 24);
  ObservationFixture f;
  cosmo_prepared *observation = nullptr;
  uint32_t semantic = 0;
  CHECK(cosmo_prepare_observations(&f.d, &f.op, &observation, &semantic) ==
            COSMO_OK &&
        observation);
  cosmo_supernova_policy policy{};
  policy.struct_size = sizeof(policy);
  policy.abi_version = COSMO_ABI_VERSION;
  policy.arithmetic = 1;
  policy.include_residual_arrays = 1;
  policy.maximum_models = 3;
  policy.maximum_source_rows = 2;
  policy.maximum_matrix_elements = 4;
  policy.maximum_array_elements = 64;
  policy.maximum_native_output_bytes = 65536;
  policy.maximum_total_evaluations = 100000;
  policy.maximum_evaluations_per_integral = 10000;
  policy.maximum_depth = 24;
  policy.absolute_tolerance = policy.relative_tolerance = 1e-12;
  policy.maximum_forward_sensitivity = 1e-10;
  cosmo_supernova *owner = nullptr;
  CHECK(cosmo_supernova_prepare(observation, &policy, &owner) == COSMO_OK &&
        owner);
  cosmo_observation_destroy(observation);
  irred::observations::Input native_input{};
  native_input.profile = Profile::pantheon_plus_released_v1;
  native_input.role = Role::released_fitted_summary;
  native_input.unit = Unit::magnitude;
  native_input.calibration = Calibration::unknown;
  native_input.uncertainty = Uncertainty::covariance;
  native_input.uncertainty_unit = UncertaintyUnit::magnitude_squared;
  native_input.table_sha256 = f.hash;
  native_input.uncertainty_sha256 = f.covhash;
  native_input.ordering_provenance = f.order;
  native_input.measurement_ids = {f.a, f.b};
  native_input.event_ids = {f.event, f.event};
  native_input.uncertainty_axis_ids = native_input.measurement_ids;
  native_input.values = {2, -3};
  native_input.zhd = {.1, .2};
  native_input.zcmb = {.1, .2};
  native_input.zhel = {.1, .2};
  native_input.missing = {0, 0};
  native_input.zhd_missing = {0, 0};
  native_input.zcmb_missing = {0, 0};
  native_input.zhel_missing = {0, 0};
  native_input.quality = {0, 0};
  native_input.uncertainty_matrix = {4, 1, 1, 9};
  auto native_source =
      irred::observations::prepare(std::move(native_input), {2, 4, 4096});
  irred::supernova::Policy np{};
  np.arithmetic = irred::numerics::Arithmetic::longdouble_cpu_v1;
  np.maximum_models = 4;
  np.maximum_matrix_elements = 4;
  np.maximum_forward_sensitivity = 1e-10;
  np.background.maximum_queries = 2;
  np.background.maximum_total_evaluations = 100000;
  np.background.integration = {1e-12, 1e-12, 10000, 24};
  auto direct = irred::supernova::prepare(std::move(native_source), np);
  CHECK(direct.status() == irred::supernova::Status::ok);
  cosmo_supernova_model_v2 parameters[] = {
      {2, 0, .3, 0, -.9, .4}, {0, 0, .3, 0, -1, 0}, {2, 0, .3, 1, -.9, .4}};
  cosmo_supernova_batch_v2 batch{sizeof(batch), COSMO_ABI_VERSION, parameters,
                                 3, sizeof(parameters)};
  cosmo_supernova_result_v2 *result = nullptr;
  CHECK(cosmo_supernova_evaluate_v2(owner, &batch, &policy, &result) ==
            COSMO_OK &&
        result);
  const cosmo_supernova_slot_v2 *rows = nullptr;
  uint64_t count = 0;
  uint32_t status = 0;
  CHECK(cosmo_supernova_result_v2_view(result, &rows, &count, &status) ==
            COSMO_OK &&
        count == 3 && status == 0);
  CHECK(rows[0].status == 0 && rows[1].status == 0 &&
        rows[0].source_parameters.w0 == -.9 &&
        rows[0].source_parameters.wa == .4);
  std::array<irred::supernova::ModelPointV2, 3> native_points{
      {{irred::cosmology::Model::flat_cpl_late_v1, .3, 0, -.9, .4},
       {irred::cosmology::Model::flat_lcdm_late_v1, .3, 0, -1, 0},
       {irred::cosmology::Model::flat_cpl_late_v1, .3, 1, -.9, .4}}};
  auto native = direct.evaluate_batch_v2(native_points, np);
  for (std::size_t i = 0; i < 3; ++i) {
    const auto &expected = native.slots[i].calculation;
    CHECK(rows[i].status == static_cast<uint32_t>(expected.status) &&
          rows[i].background_status ==
              static_cast<uint32_t>(expected.background_status) &&
          rows[i].numerical_status ==
              static_cast<uint32_t>(expected.numerical_status));
    if (expected.status == irred::supernova::Status::ok) {
      CHECK(rows[i].quadratic == expected.quadratic &&
            rows[i].relative_profile_score == expected.relative_profile_score &&
            rows[i].offset_coefficient == expected.offset_coefficient);
      CHECK(rows[i].shape_magnitudes.length ==
            expected.shape_magnitudes.size());
      for (std::size_t j = 0; j < expected.shape_magnitudes.size(); ++j)
        CHECK(rows[i].shape_magnitudes.data[j] ==
                  expected.shape_magnitudes[j] &&
              rows[i].profiled_residuals.data[j] ==
                  expected.profiled_residuals[j]);
    }
  }
  CHECK(rows[2].status != 0 && rows[2].source_parameters.constant_q == 1 &&
        rows[2].source_parameters.wa == .4 &&
        rows[2].shape_magnitudes.length == 0 &&
        rows[2].profiled_residuals.length == 0);
  cosmo_supernova_model legacy{0, 0, .3, 0};
  cosmo_supernova_batch old_batch{sizeof(old_batch), COSMO_ABI_VERSION, &legacy,
                                  1, sizeof(legacy)};
  cosmo_supernova_result *old_result = nullptr;
  CHECK(cosmo_supernova_evaluate(owner, &old_batch, &policy, &old_result) ==
            COSMO_OK &&
        old_result);
  const cosmo_supernova_slot *oldrows = nullptr;
  uint64_t oldcount = 0;
  uint32_t oldstatus = 0;
  CHECK(cosmo_supernova_result_view(old_result, &oldrows, &oldcount,
                                    &oldstatus) == COSMO_OK &&
        oldrows[0].relative_profile_score == rows[1].relative_profile_score &&
        oldrows[0].shape_magnitudes.data[0] ==
            rows[1].shape_magnitudes.data[0]);
  cosmo_supernova_result_destroy(old_result);
  parameters[0].w0 = -2;
  CHECK(rows[0].source_parameters.w0 == -.9);
  parameters[0].w0 = -.9;
  cosmo_supernova_result_v2_destroy(result);
  auto before = live;
  int failures = 0;
  bool reached = false;
  for (int k = 0; k < 200; ++k) {
    result = nullptr;
    fail_after = k;
    auto code = cosmo_supernova_evaluate_v2(owner, &batch, &policy, &result);
    fail_after = -1;
    if (code == COSMO_OK) {
      CHECK(result);
      cosmo_supernova_result_v2_destroy(result);
      CHECK(live == before);
      reached = true;
      break;
    }
    CHECK(code == COSMO_ALLOCATION_FAILURE && !result && live == before);
    ++failures;
  }
  CHECK(reached && failures > 0);
  std::printf("supernova_v2 allocation_failures %d\n", failures);
  auto impossible = batch;
  impossible.models = nullptr;
  impossible.model_count = UINT64_MAX;
  impossible.model_byte_length = 0;
  fail_after = 0;
  auto code = cosmo_supernova_evaluate_v2(owner, &impossible, &policy, &result);
  fail_after = -1;
  CHECK(code == COSMO_INVALID_INPUT && !result && live == before);
  auto wrong = batch;
  wrong.model_byte_length--;
  CHECK(cosmo_supernova_evaluate_v2(owner, &wrong, &policy, &result) ==
            COSMO_INVALID_INPUT &&
        !result);
  wrong = batch;
  wrong.abi_version++;
  CHECK(cosmo_supernova_evaluate_v2(owner, &wrong, &policy, &result) ==
            COSMO_ABI_MISMATCH &&
        !result);
  parameters[0].reserved = 1;
  CHECK(cosmo_supernova_evaluate_v2(owner, &batch, &policy, &result) ==
            COSMO_INVALID_INPUT &&
        !result);
  parameters[0].reserved = 0;
  parameters[1].wa = .1;
  CHECK(cosmo_supernova_evaluate_v2(owner, &batch, &policy, &result) ==
            COSMO_OK &&
        result);
  CHECK(cosmo_supernova_result_v2_view(result, &rows, &count, &status) ==
            COSMO_OK &&
        rows[1].status != 0 && rows[1].source_parameters.wa == .1 &&
        rows[1].profiled_residuals.length == 0);
  cosmo_supernova_result_v2_destroy(result);
  parameters[1].wa = 0;
  CHECK(cosmo_supernova_evaluate_v2(owner, &batch, &policy, nullptr) ==
        COSMO_INVALID_INPUT);
  wrong = batch;
  wrong.models = reinterpret_cast<const cosmo_supernova_model_v2 *>(
      reinterpret_cast<const char *>(parameters) + 1);
  CHECK(cosmo_supernova_evaluate_v2(owner, &wrong, &policy, &result) ==
            COSMO_INVALID_INPUT &&
        !result);
  auto cap = policy;
  cap.maximum_array_elements = 1;
  CHECK(cosmo_supernova_evaluate_v2(owner, &batch, &cap, &result) ==
            COSMO_INVALID_INPUT &&
        !result);
  cap = policy;
  cap.maximum_total_evaluations = 3;
  CHECK(cosmo_supernova_evaluate_v2(owner, &batch, &cap, &result) == COSMO_OK &&
        result);
  CHECK(cosmo_supernova_result_v2_view(result, &rows, &count, &status) ==
        COSMO_OK);
  uint64_t calls = 0;
  for (uint64_t i = 0; i < count; ++i)
    calls += rows[i].background_evaluations;
  CHECK(calls <= 3 && rows[0].status != 0 && rows[1].status != 0);
  cosmo_supernova_result_v2_destroy(result);
  cap = policy;
  cap.maximum_forward_sensitivity = 1e-30;
  CHECK(cosmo_supernova_evaluate_v2(owner, &batch, &cap, &result) == COSMO_OK &&
        result);
  CHECK(cosmo_supernova_result_v2_view(result, &rows, &count, &status) ==
            COSMO_OK &&
        rows[0].numerical_status ==
            COSMO_NUMERICAL_STATUS_CONDITIONING_BUDGET_EXCEEDED &&
        rows[0].shape_magnitudes.length == 0);
  cosmo_supernova_result_v2_destroy(result);
  CHECK(cosmo_supernova_evaluate_v2(owner, &batch, &policy, &result) ==
            COSMO_OK &&
        result);
  cosmo_supernova_destroy(owner);
  CHECK(cosmo_supernova_result_v2_view(result, &rows, &count, &status) ==
            COSMO_OK &&
        rows[0].shape_magnitudes.data && rows[0].source_parameters.wa == .4);
  cosmo_supernova_result_v2_destroy(result);
}
int main() {
  try {
    background_test();
    supernova_test();
    CHECK(live == 0);
    std::printf("CPL v2 ABI %d PASS\n", checks);
    return 0;
  } catch (const std::exception &e) {
    fail_after = -1;
    std::fprintf(stderr, "FAIL %s checks%d\n", e.what(), checks);
    return 1;
  }
}
