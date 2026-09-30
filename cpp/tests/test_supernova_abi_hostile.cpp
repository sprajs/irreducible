// Transport adversaries by consumer author, independently of bridge
// implementation. Generated release-shaped bytes test interface semantics, not
// actual survey science.
#include "irred/abi.h"
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
struct Fixture {
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
  Fixture() {
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
int main() {
  try {
    Fixture f;
    cosmo_prepared *obs = nullptr;
    uint32_t semantic = 99;
    CHECK(cosmo_prepare_observations(&f.d, &f.op, &obs, &semantic) ==
              COSMO_OK &&
          obs);
    cosmo_supernova_policy p{};
    p.struct_size = sizeof(p);
    p.abi_version = COSMO_ABI_VERSION;
    p.arithmetic = 1;
    p.include_residual_arrays = 1;
    p.maximum_models = 4;
    p.maximum_source_rows = 2;
    p.maximum_matrix_elements = 4;
    p.maximum_array_elements = 64;
    p.maximum_native_output_bytes = 65536;
    p.maximum_total_evaluations = 100000;
    p.maximum_evaluations_per_integral = 10000;
    p.maximum_depth = 24;
    p.absolute_tolerance = 1e-12;
    p.relative_tolerance = 1e-12;
    p.maximum_forward_sensitivity = 1e-10;
    cosmo_supernova *owner = nullptr;
    CHECK(cosmo_supernova_prepare(obs, &p, &owner) == COSMO_OK && owner);
    cosmo_supernova_view view{};
    CHECK(cosmo_supernova_source_view(owner, &view) == COSMO_OK &&
          view.status == 0 && view.ordered_ids.length == 2);
    CHECK(view.selected_z_expansion.data[0] == .1 &&
          view.selected_z_observer.data[1] == .2);
    CHECK(view.selected_source_indices.length == 2 &&
          view.selected_source_indices.data[0] == 0 &&
          view.selected_source_indices.data[1] == 1);
    CHECK(std::string(
              reinterpret_cast<const char *>(view.ordered_ids.data[0].data),
              view.ordered_ids.data[0].length) == "a");
    CHECK(std::string(
              reinterpret_cast<const char *>(view.ordered_ids.data[1].data),
              view.ordered_ids.data[1].length) == "b");
    bool reached = false;
    int failures = 0;
    for (int k = 0; k < 300; ++k) {
      auto before = live;
      cosmo_supernova *tmp = nullptr;
      fail_after = k;
      auto st = cosmo_supernova_prepare(obs, &p, &tmp);
      fail_after = -1;
      if (st == COSMO_OK) {
        CHECK(tmp);
        cosmo_supernova_destroy(tmp);
        CHECK(live == before);
        reached = true;
        break;
      }
      CHECK(st == COSMO_ALLOCATION_FAILURE && !tmp && live == before);
      ++failures;
    }
    CHECK(reached && failures > 0);
    const int prepare_failures = failures;
    auto bad = p;
    bad.arithmetic = 99;
    cosmo_supernova *tmp = nullptr;
    CHECK(cosmo_supernova_prepare(obs, &bad, &tmp) == COSMO_INVALID_INPUT &&
          !tmp);
    bad = p;
    bad.abi_version++;
    CHECK(cosmo_supernova_prepare(obs, &bad, &tmp) == COSMO_ABI_MISMATCH &&
          !tmp);
    bad = p;
    bad.maximum_source_rows = 1;
    CHECK(cosmo_supernova_prepare(obs, &bad, &tmp) == COSMO_INVALID_INPUT &&
          !tmp);
    bad = p;
    bad.maximum_forward_sensitivity = 1e-30;
    cosmo_supernova *failed_owner = nullptr;
    CHECK(cosmo_supernova_prepare(obs, &bad, &failed_owner) == COSMO_OK &&
          failed_owner);
    cosmo_supernova_view failed_view{};
    CHECK(cosmo_supernova_source_view(failed_owner, &failed_view) == COSMO_OK &&
          failed_view.status != 0);
    CHECK(failed_view.preparation_numerical_status ==
              COSMO_NUMERICAL_STATUS_CONDITIONING_BUDGET_EXCEEDED &&
          failed_view.matrix_validation_assessed == 0 &&
          failed_view.arithmetic_id.length == 0);
    CHECK(failed_view.source_row_count == 2 &&
          failed_view.table_identity.length == 64);
    CHECK(cosmo_supernova_destroy(failed_owner) == COSMO_OK);
    CHECK(cosmo_observation_destroy(obs) == COSMO_OK);
    obs = nullptr;
    f.values[0] = 999;
    f.z[0] = 4;
    CHECK(cosmo_supernova_source_view(owner, &view) == COSMO_OK &&
          view.selected_z_expansion.data[0] == .1);
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
    std::array<irred::supernova::ModelPoint, 2> native_points{
        {{irred::cosmology::Model::flat_lcdm_late_v1, .3, 0},
         {irred::cosmology::Model::constant_q_flat_v1, 0, 0}}};
    auto direct_batch = direct.evaluate_batch(native_points, np);
    cosmo_supernova_model points[2] = {{0, 0, .3, 0}, {1, 0, 0, 0}};
    cosmo_supernova_batch batch{sizeof(batch), COSMO_ABI_VERSION, points, 2,
                                sizeof(points)};
    cosmo_supernova_result *result = nullptr;
    CHECK(cosmo_supernova_evaluate(owner, &batch, &p, &result) == COSMO_OK &&
          result);
    const cosmo_supernova_slot *rows = nullptr;
    uint64_t count = 0;
    uint32_t status = 99;
    CHECK(cosmo_supernova_result_view(result, &rows, &count, &status) ==
              COSMO_OK &&
          count == 2 && status == 0);
    CHECK(rows[0].model_index == 0 && rows[1].model_index == 1 &&
          rows[0].status == 0 && rows[1].status == 0);
    CHECK(rows[0].shape_magnitudes.length == 2 &&
          rows[0].profiled_residuals.length == 2);
    CHECK(rows[0].relative_profile_score ==
              direct_batch.slots[0].relative_profile_score &&
          rows[1].quadratic == direct_batch.slots[1].quadratic);
    CHECK(rows[0].offset_coefficient ==
          direct_batch.slots[0].offset_coefficient);
    CHECK(rows[0].profiled_residuals.data[0] ==
          direct_batch.slots[0].profiled_residuals[0]);
    CHECK(std::abs(rows[0].relative_profile_score + .5 * rows[0].quadratic) <
          1e-12);
    CHECK(cosmo_supernova_result_destroy(result) == COSMO_OK);
    reached = false;
    failures = 0;
    for (int k = 0; k < 300; ++k) {
      auto before = live;
      result = nullptr;
      fail_after = k;
      auto st = cosmo_supernova_evaluate(owner, &batch, &p, &result);
      fail_after = -1;
      if (st == COSMO_OK) {
        CHECK(result);
        cosmo_supernova_result_destroy(result);
        CHECK(live == before);
        reached = true;
        break;
      }
      CHECK(st == COSMO_ALLOCATION_FAILURE && !result && live == before);
      ++failures;
    }
    CHECK(reached && failures > 0);
    auto impossible = batch;
    impossible.models = nullptr;
    impossible.model_count = UINT64_MAX;
    impossible.model_byte_length = 0;
    auto before_invalid = live;
    fail_after = 0;
    auto rejected = cosmo_supernova_evaluate(owner, &impossible, &p, &result);
    fail_after = -1;
    CHECK(rejected == COSMO_INVALID_INPUT && !result && live == before_invalid);
    auto misaligned = batch;
    misaligned.models = reinterpret_cast<const cosmo_supernova_model *>(
        reinterpret_cast<const char *>(points) + 1);
    CHECK(cosmo_supernova_evaluate(owner, &misaligned, &p, &result) ==
              COSMO_INVALID_INPUT &&
          !result);
    CHECK(cosmo_supernova_evaluate(nullptr, &batch, &p, &result) ==
              COSMO_INVALID_INPUT &&
          !result);
    auto malformed = batch;
    malformed.model_byte_length--;
    CHECK(cosmo_supernova_evaluate(owner, &malformed, &p, &result) ==
              COSMO_INVALID_INPUT &&
          !result);
    bad = p;
    bad.maximum_array_elements = 1;
    CHECK(cosmo_supernova_evaluate(owner, &batch, &bad, &result) ==
              COSMO_INVALID_INPUT &&
          !result);
    bad = p;
    bad.maximum_native_output_bytes = 1;
    CHECK(cosmo_supernova_evaluate(owner, &batch, &bad, &result) ==
              COSMO_INVALID_INPUT &&
          !result);
    bad = p;
    bad.maximum_total_evaluations = 3;
    CHECK(cosmo_supernova_evaluate(owner, &batch, &bad, &result) == COSMO_OK &&
          result);
    CHECK(cosmo_supernova_result_view(result, &rows, &count, &status) ==
              COSMO_OK &&
          count == 2);
    CHECK(rows[0].background_evaluations + rows[1].background_evaluations <= 3);
    CHECK(rows[0].status != 0 && rows[1].status != 0 &&
          rows[1].shape_magnitudes.length == 0);
    cosmo_supernova_result_destroy(result);
    points[1].model = 99;
    CHECK(cosmo_supernova_evaluate(owner, &batch, &p, &result) == COSMO_OK &&
          result);
    CHECK(cosmo_supernova_result_view(result, &rows, &count, &status) ==
              COSMO_OK &&
          count == 2);
    CHECK(rows[1].status != 0 && rows[1].shape_magnitudes.length == 0 &&
          rows[1].profiled_residuals.length == 0);
    CHECK(cosmo_supernova_destroy(owner) == COSMO_OK);
    owner = nullptr;
    CHECK(cosmo_supernova_result_view(result, &rows, &count, &status) ==
              COSMO_OK &&
          count == 2 && rows[0].shape_magnitudes.length == 2);
    cosmo_supernova_result_destroy(result);
    CHECK(cosmo_supernova_destroy(owner) == COSMO_OK);
    std::printf("supernova ABI hostile: %d checks PASS; prepare allocation "
                "failures=%d, evaluate allocation failures=%d\n",
                checks, prepare_failures, failures);
  } catch (const std::exception &e) {
    fail_after = -1;
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
