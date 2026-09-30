// Independent BAO transport adversaries; wrapper parity is not independent
// science.
#include "irred/abi.h"
#include "irred/background.hpp"
#include "irred/bao.hpp"
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
      throw std::runtime_error("BAO ABI: " #__VA_ARGS__);                      \
  } while (false)
using namespace irred;
cosmo_bytes bytes(const std::string &s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
cosmo_f64_buffer doubles(const double *p, std::size_t n) {
  return {sizeof(cosmo_f64_buffer), COSMO_ABI_VERSION, 2, 0, p, n,
          n * sizeof(double)};
}
struct Fixture {
  std::string a = "synthetic:DM:.1", b = "synthetic:DH:.2",
              table = "inline synthetic table",
              cov = "rational covariance det35",
              order = "supplied exact two-coordinate order",
              cal = "synthetic no calibration",
              dep = "correlated control; cross-probe dependence unknown",
              red = "P01/released-effective-redshift/v1",
              ruler = "P01/free-H0rd-km-s-no-early-physics/v1",
              h0 = "P01/computational-H0-fixed-70-km-s-Mpc/v1";
  cosmo_bytes ids[2] = {bytes(a), bytes(b)};
  cosmo_bao_query queries[2] = {{.1, 0, 0}, {.2, 1, 0}};
  double observed[2] = {299792.458 / 10000 * .1 + 1, 299792.458 / 10000 + 2};
  double matrix[4] = {4, 1, 1, 9};
  cosmo_bao_descriptor descriptor{};
  cosmo_bao_policy policy{};
  Fixture() {
    descriptor.struct_size = sizeof(descriptor);
    descriptor.abi_version = COSMO_ABI_VERSION;
    descriptor.role = COSMO_BAO_ROLE_SYNTHETIC_CONTROL;
    descriptor.covariance_unit = COSMO_BAO_COVARIANCE_UNIT_RATIO_SQUARED;
    descriptor.queries = queries;
    descriptor.query_count = 2;
    descriptor.query_byte_length = sizeof queries;
    descriptor.observed = doubles(observed, 2);
    descriptor.covariance = doubles(matrix, 4);
    descriptor.ordered_ids =
        descriptor.covariance_axis_ids = {ids, 2, sizeof ids};
    descriptor.table_identity = bytes(table);
    descriptor.covariance_identity = bytes(cov);
    descriptor.ordering_provenance = bytes(order);
    descriptor.calibration_provenance = bytes(cal);
    descriptor.dependence_provenance = bytes(dep);
    descriptor.redshift_convention = bytes(red);
    descriptor.ruler_convention = bytes(ruler);
    descriptor.computational_h0_convention = bytes(h0);
    policy.struct_size = sizeof(policy);
    policy.abi_version = COSMO_ABI_VERSION;
    policy.background = {sizeof(cosmo_background_policy),
                         COSMO_ABI_VERSION,
                         30,
                         0,
                         8,
                         2,
                         16,
                         1000000,
                         2000000,
                         100000,
                         1e-13,
                         1e-12};
    policy.arithmetic = 1;
    policy.include_predictions = policy.include_residuals = 1;
    policy.maximum_models = 8;
    policy.maximum_rows = 2;
    policy.maximum_matrix_elements = 4;
    policy.maximum_string_bytes = 4096;
    policy.maximum_native_bytes = 1000000;
    policy.maximum_array_elements = 32;
    policy.maximum_native_output_bytes = 1000000;
    policy.maximum_forward_sensitivity = 1e-10;
  }
  bao::DensityInput native_input() const {
    bao::DensityInput n;
    n.queries = {{.1, bao::Observable::transverse_over_ruler},
                 {.2, bao::Observable::hubble_over_ruler}};
    n.observed = {observed[0], observed[1]};
    n.covariance = {4, 1, 1, 9};
    n.ordered_ids = {a, b};
    n.role = bao::RowRole::synthetic_control;
    n.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
    n.table_identity = table;
    n.covariance_identity = cov;
    n.ordering_provenance = order;
    n.calibration_provenance = cal;
    n.dependence_provenance = dep;
    return n;
  }
  bao::DensityPolicy native_policy() const {
    bao::Policy o;
    o.maximum_queries = 2;
    o.maximum_native_bytes = 1000000;
    o.background = {{1e-13, 1e-12, 100000, 30}, 2, 2000000};
    return {
        o, 8, 4, 4096, 1000000, 1e-10, numerics::Arithmetic::longdouble_cpu_v1};
  }
};
std::string string(cosmo_bytes b) {
  return {reinterpret_cast<const char *>(b.data),
          static_cast<std::size_t>(b.length)};
}
int main() {
  try {
    Fixture f;
    static_assert(sizeof(cosmo_bao_query) == 16);
    static_assert(offsetof(cosmo_bao_model, h0_rd_km_s) == 40);
    cosmo_bao *owner = nullptr;
    CHECK(cosmo_bao_prepare(&f.descriptor, &f.policy, &owner) == COSMO_OK &&
          owner);
    cosmo_bao_source_view_t view{};
    view.struct_size = sizeof(view);
    view.abi_version = COSMO_ABI_VERSION;
    CHECK(cosmo_bao_source_view(owner, &view) == COSMO_OK && view.status == 0);
    CHECK(view.source.query_count == 2 && view.source.role == 1 &&
          view.source.covariance.data[1] == 1);
    CHECK(string(view.source.ordering_provenance) == f.order &&
          string(view.source.redshift_convention) == f.red);
    CHECK(view.prepare_policy.maximum_forward_sensitivity == 1e-10 &&
          view.prepare_policy.arithmetic == 1);
    cosmo_bao_model models[] = {{{1, 0, 0, -1, -1, 0}, 10000},
                                {{0, 0, .3, 0, -1, 0}, 5000},
                                {{2, 0, .3, 0, -.9, .4}, 15000},
                                {{2, 0, .3, 1, -.9, .4}, 10000}};
    cosmo_bao_batch batch{sizeof(batch), COSMO_ABI_VERSION, models, 4,
                          sizeof models};
    cosmo_bao_result *result = nullptr;
    CHECK(cosmo_bao_evaluate(owner, &batch, &f.policy, &result) == COSMO_OK &&
          result);
    const cosmo_bao_row *rows = nullptr;
    uint64_t count = 0, work = 0;
    uint32_t status = 99, numerical = 99;
    CHECK(cosmo_bao_result_view(result, &rows, &count, &status, &numerical,
                                &work) == COSMO_OK &&
          count == 4);
    auto native = bao::prepare_density(f.native_input(), f.native_policy());
    std::vector<bao::ModelQuery> points{
        {cosmology::prepare({cosmology::Model::constant_q_flat_v1, 70, 0, -1}),
         bao::Ruler(10000)},
        {cosmology::prepare({cosmology::Model::flat_lcdm_late_v1, 70, .3, 0}),
         bao::Ruler(5000)},
        {cosmology::prepare_cpl({70, .3, -.9, .4}), bao::Ruler(15000)}};
    auto expected = native.evaluate(points, f.native_policy());
    for (std::size_t i = 0; i < 3; ++i) {
      auto &r = rows[i];
      auto &n = expected.slots[i];
      CHECK(r.model_index == i && r.status == 0 && r.numerical_status == 0);
      CHECK(r.log_density == n.result.density.log_value &&
            r.quadratic == n.result.quadratic &&
            r.log_determinant == n.result.log_determinant &&
            r.normalization == n.result.normalization);
      CHECK(r.predictions.length == 2 && r.residuals.length == 2);
      for (std::size_t j = 0; j < 2; ++j)
        CHECK(r.predictions.data[j] == n.predictions[j] &&
              r.residuals.data[j] == n.residuals[j]);
      CHECK(r.source_parameters.h0_rd_km_s == models[i].h0_rd_km_s);
    }
    CHECK(std::abs(rows[0].quadratic - .6) < 1e-13 &&
          std::abs(rows[0].log_density +
                   .5 * (.6 + std::log(35.) +
                         2 * std::log(2 * std::acos(-1.)))) < 1e-13);
    CHECK(rows[3].status != 0 &&
          rows[3].source_parameters.parameters.constant_q == 1 &&
          rows[3].predictions.length == 0 && rows[3].residuals.length == 0);
    f.observed[0] = 999;
    f.matrix[1] = 999;
    f.queries[0].z = 5;
    models[0].h0_rd_km_s = 0;
    CHECK(rows[0].source_parameters.h0_rd_km_s == 10000 &&
          view.source.observed.data[0] != 999 &&
          view.source.queries[0].z == .1);
    CHECK(cosmo_bao_destroy(owner) == COSMO_OK);
    owner = nullptr;
    CHECK(rows[0].ordered_ids.length == 2 &&
          string(rows[0].ordered_ids.data[0]) == f.a &&
          rows[0].query_count == 2 && rows[0].queries[0].z == .1);
    CHECK(rows[0].predictions.data[0] == expected.slots[0].predictions[0] &&
          string(rows[0].ruler_convention_id) == f.ruler);
    CHECK(cosmo_bao_result_destroy(result) == COSMO_OK);
    result = nullptr;
    f.observed[0] = 299792.458 / 10000 * .1 + 1;
    f.matrix[1] = 1;
    f.queries[0].z = .1;
    models[0].h0_rd_km_s = 10000;
    // Guard descriptors before reading oversized caller ranges.
    auto d = f.descriptor;
    d.query_count = UINT64_MAX;
    d.query_byte_length = UINT64_MAX;
    auto base = live;
    CHECK(cosmo_bao_prepare(&d, &f.policy, &owner) != COSMO_OK && !owner &&
          live == base);
    d = f.descriptor;
    d.query_byte_length--;
    CHECK(cosmo_bao_prepare(&d, &f.policy, &owner) != COSMO_OK && !owner);
    d = f.descriptor;
    d.abi_version++;
    CHECK(cosmo_bao_prepare(&d, &f.policy, &owner) == COSMO_ABI_MISMATCH &&
          !owner);
    d = f.descriptor;
    d.observed.data = reinterpret_cast<const double *>(
        reinterpret_cast<const char *>(f.observed) + 1);
    CHECK(cosmo_bao_prepare(&d, &f.policy, &owner) != COSMO_OK && !owner);
    // Complete allocation sweeps, stable errors and cleanup in original
    // allocator.
    for (int field = 0; field < 3; ++field) {
      d = f.descriptor;
      if (field == 0)
        d.role = UINT32_MAX;
      if (field == 1)
        d.covariance_unit = UINT32_MAX;
      if (field == 2)
        d.redshift_convention = bytes(f.dep);
      base = live;
      CHECK(cosmo_bao_prepare(&d, &f.policy, &owner) == COSMO_INVALID_INPUT &&
            !owner && live == base);
    }
    cosmo_bytes reordered[2] = {f.ids[1], f.ids[0]};
    d = f.descriptor;
    d.covariance_axis_ids = {reordered, 2, sizeof reordered};
    CHECK(cosmo_bao_prepare(&d, &f.policy, &owner) == COSMO_INVALID_INPUT &&
          !owner);
    f.queries[0].observable = UINT32_MAX;
    CHECK(cosmo_bao_prepare(&f.descriptor, &f.policy, &owner) == COSMO_OK &&
          owner);
    CHECK(cosmo_bao_source_view(owner, &view) == COSMO_OK &&
          view.status != COSMO_GAUSSIAN_STATUS_FINITE &&
          view.source.queries[0].observable == UINT32_MAX);
    CHECK(cosmo_bao_destroy(owner) == COSMO_OK);
    owner = nullptr;
    f.queries[0].observable = 0;
    f.matrix[1] = 2;
    CHECK(cosmo_bao_prepare(&f.descriptor, &f.policy, &owner) == COSMO_OK &&
          owner);
    CHECK(cosmo_bao_source_view(owner, &view) == COSMO_OK &&
          view.status != COSMO_GAUSSIAN_STATUS_FINITE &&
          view.source.covariance.data[1] == 2 &&
          view.source.covariance.data[2] == 1);
    CHECK(cosmo_bao_evaluate(owner, &batch, &f.policy, &result) == COSMO_OK &&
          result);
    CHECK(cosmo_bao_result_view(result, &rows, &count, &status, &numerical,
                                &work) == COSMO_OK &&
          count == 0 && status != COSMO_GAUSSIAN_STATUS_FINITE);
    CHECK(cosmo_bao_result_destroy(result) == COSMO_OK);
    result = nullptr;
    CHECK(cosmo_bao_destroy(owner) == COSMO_OK);
    owner = nullptr;
    f.matrix[1] = 1;
    int prepare_failures = 0, evaluate_failures = 0;
    bool reached = false;
    for (int i = 0; i < 512; ++i) {
      base = live;
      fail_after = i;
      auto call = cosmo_bao_prepare(&f.descriptor, &f.policy, &owner);
      fail_after = -1;
      if (call == COSMO_OK) {
        CHECK(owner);
        CHECK(cosmo_bao_destroy(owner) == COSMO_OK);
        owner = nullptr;
        CHECK(live == base);
        reached = true;
        break;
      }
      CHECK(call == COSMO_ALLOCATION_FAILURE && !owner && live == base);
      ++prepare_failures;
    }
    CHECK(reached && prepare_failures > 0);
    CHECK(cosmo_bao_prepare(&f.descriptor, &f.policy, &owner) == COSMO_OK &&
          owner);
    auto no_arrays = f.policy;
    no_arrays.include_predictions = no_arrays.include_residuals = 0;
    CHECK(cosmo_bao_evaluate(owner, &batch, &no_arrays, &result) == COSMO_OK &&
          result);
    CHECK(cosmo_bao_result_view(result, &rows, &count, &status, &numerical,
                                &work) == COSMO_OK &&
          count == 4);
    CHECK(rows[0].status == 0 &&
          rows[0].log_density == expected.slots[0].result.density.log_value &&
          rows[0].predictions.length == 0 && rows[0].residuals.length == 0);
    CHECK(cosmo_bao_result_destroy(result) == COSMO_OK);
    result = nullptr;
    no_arrays.maximum_array_elements = 15;
    base = live;
    CHECK(cosmo_bao_evaluate(owner, &batch, &no_arrays, &result) != COSMO_OK &&
          !result && live == base);
    auto cap = f.policy;
    cap.background.maximum_total_evaluations = expected.slots[0].evaluations;
    CHECK(cosmo_bao_evaluate(owner, &batch, &cap, &result) == COSMO_OK &&
          result);
    CHECK(cosmo_bao_result_view(result, &rows, &count, &status, &numerical,
                                &work) == COSMO_OK &&
          count == 4);
    CHECK(rows[0].status == 0 && rows[1].status != 0 &&
          rows[1].numerical_status == COSMO_NUMERICAL_STATUS_WORK_LIMIT &&
          work <= cap.background.maximum_total_evaluations);
    CHECK(rows[1].predictions.length == 0 && rows[1].residuals.length == 0 &&
          rows[1].log_density == 0);
    CHECK(cosmo_bao_result_destroy(result) == COSMO_OK);
    result = nullptr;
    auto tight = f.policy;
    tight.maximum_forward_sensitivity = 1e-30;
    CHECK(cosmo_bao_evaluate(owner, &batch, &tight, &result) == COSMO_OK &&
          result);
    CHECK(cosmo_bao_result_view(result, &rows, &count, &status, &numerical,
                                &work) == COSMO_OK);
    CHECK(rows[0].status != 0 &&
          rows[0].numerical_status ==
              COSMO_NUMERICAL_STATUS_CONDITIONING_BUDGET_EXCEEDED &&
          rows[0].predictions.length == 0);
    CHECK(cosmo_bao_result_destroy(result) == COSMO_OK);
    result = nullptr;
    auto bad_batch = batch;
    bad_batch.model_count = UINT64_MAX;
    bad_batch.model_byte_length = UINT64_MAX;
    base = live;
    CHECK(cosmo_bao_evaluate(owner, &bad_batch, &f.policy, &result) !=
              COSMO_OK &&
          !result && live == base);
    bad_batch = batch;
    bad_batch.model_byte_length--;
    CHECK(cosmo_bao_evaluate(owner, &bad_batch, &f.policy, &result) !=
              COSMO_OK &&
          !result);
    cap = f.policy;
    cap.maximum_native_output_bytes = 1;
    base = live;
    CHECK(cosmo_bao_evaluate(owner, &batch, &cap, &result) != COSMO_OK &&
          !result && live == base);
    reached = false;
    for (int i = 0; i < 512; ++i) {
      base = live;
      fail_after = i;
      auto call = cosmo_bao_evaluate(owner, &batch, &f.policy, &result);
      fail_after = -1;
      if (call == COSMO_OK) {
        CHECK(result);
        CHECK(cosmo_bao_result_destroy(result) == COSMO_OK);
        result = nullptr;
        CHECK(live == base);
        reached = true;
        break;
      }
      CHECK(call == COSMO_ALLOCATION_FAILURE && !result && live == base);
      ++evaluate_failures;
    }
    CHECK(reached && evaluate_failures > 0);
    CHECK(cosmo_bao_destroy(owner) == COSMO_OK);
    owner = nullptr;
    std::printf("BAO ABI checks %d PASS allocation_sweep %d/%d\n", checks,
                prepare_failures, evaluate_failures);
    return 0;
  } catch (const std::exception &e) {
    fail_after = -1;
    std::fprintf(stderr, "FAIL %s after %d\n", e.what(), checks);
    return 1;
  }
}
