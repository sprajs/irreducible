// Independent BAO transport adversaries; wrapper parity is not independent
// science.
#include "irred/abi.h"
#include "irred/background.hpp"
#include "irred/bao.hpp"
#include "irred/supernova.hpp"
#include <algorithm>
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
  std::string a = "synthetic:DM:.1", b = std::string(1024, 'L'),
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
cosmo_bao_piecewise_policy analytic(const Fixture &f) {
  cosmo_bao_piecewise_policy p{};
  p.struct_size = sizeof(p);
  p.abi_version = COSMO_ABI_VERSION;
  p.arithmetic = 1;
  p.include_predictions = p.include_residuals = 1;
  p.maximum_models = 8;
  p.maximum_rows = 2;
  p.maximum_matrix_elements = 4;
  p.maximum_string_bytes = 8192;
  p.maximum_native_bytes = p.maximum_native_output_bytes = 1000000;
  p.maximum_array_elements = 32;
  p.maximum_queries = 2;
  p.maximum_total_segment_visits = 2000;
  p.maximum_forward_sensitivity = f.policy.maximum_forward_sensitivity;
  return p;
}
bao::PiecewiseDensityPolicy native_policy(const cosmo_bao_piecewise_policy &p) {
  bao::PiecewiseDensityPolicy n;
  n.arithmetic = static_cast<numerics::Arithmetic>(p.arithmetic);
  n.maximum_models = p.maximum_models;
  n.maximum_native_bytes = p.maximum_native_output_bytes;
  n.maximum_forward_sensitivity = p.maximum_forward_sensitivity;
  n.observables.maximum_queries = p.maximum_queries;
  n.observables.maximum_native_bytes = p.maximum_native_output_bytes;
  n.observables.background.maximum_queries = p.maximum_queries;
  n.observables.background.maximum_segment_visits =
      p.maximum_total_segment_visits;
  return n;
}
int main() {
  try {
    Fixture f;
    auto p = analytic(f);
    static_assert(sizeof(cosmo_bao_piecewise_model) == 48);
    static_assert(sizeof(cosmo_bao_piecewise_policy) == 104);
    static_assert(sizeof(cosmo_bao_piecewise_batch) == 32);
    static_assert(offsetof(cosmo_bao_piecewise_model, h0_rd_km_s) == 40);
    cosmo_bao_piecewise *owner = nullptr;
    CHECK(cosmo_bao_piecewise_prepare(&f.descriptor, &p, &owner) == COSMO_OK &&
          owner);
    cosmo_bao_piecewise_source_view_t source{};
    source.struct_size = sizeof(source);
    source.abi_version = COSMO_ABI_VERSION;
    CHECK(cosmo_bao_piecewise_source_view(owner, &source) == COSMO_OK &&
          source.status == 0);
    CHECK(source.source.query_count == 2 &&
          source.source.covariance.data[1] == 1 &&
          string(source.source.ordered_ids.data[1]) == f.b);
    CHECK(source.prepare_policy.maximum_queries == 2 &&
          string(source.arithmetic_id).find("longdouble") != std::string::npos);
    cosmo_bao_piecewise_model models[] = {{{-1, -1, -1, -1, -1}, 10000},
                                          {{-.4, .2, -.3, 1, -2}, 5000},
                                          {{2, 2, 2, 2, 2}, 15000},
                                          {{3, 0, 0, 0, 0}, 10000},
                                          {{0, 0, 0, 0, 0}, 0}};
    cosmo_bao_piecewise_batch batch{sizeof(batch), COSMO_ABI_VERSION, models, 5,
                                    sizeof(models)};
    cosmo_bao_piecewise_result *result = nullptr;
    CHECK(cosmo_bao_piecewise_evaluate(owner, &batch, &p, &result) ==
              COSMO_OK &&
          result);
    const cosmo_bao_piecewise_row *rows = nullptr;
    uint64_t count = 0, work = 0;
    uint32_t status = 99, cause = 99;
    CHECK(cosmo_bao_piecewise_result_view(result, &rows, &count, &status,
                                          &cause, &work) == COSMO_OK &&
          count == 5);
    auto native = bao::prepare_density(f.native_input(), f.native_policy());
    std::vector<bao::PiecewiseModelPoint> points;
    for (auto &m : models) {
      std::array<double, 5> q;
      std::copy_n(m.q, 5, q.begin());
      points.emplace_back(q, bao::Ruler(m.h0_rd_km_s));
    }
    auto expected = native.evaluate_piecewise(points, native_policy(p));
    CHECK(work == expected.segment_visits &&
          status == static_cast<uint32_t>(expected.status) &&
          cause == static_cast<uint32_t>(expected.numerical_status));
    for (std::size_t i = 0; i < count; ++i) {
      auto &r = rows[i];
      auto &e = expected.slots[i];
      CHECK(r.model_index == i &&
            r.status == static_cast<uint32_t>(e.result.density.status) &&
            r.numerical_status == static_cast<uint32_t>(e.numerical_status));
      CHECK(r.segment_visits == e.segment_visits &&
            r.source_parameters.h0_rd_km_s == models[i].h0_rd_km_s);
      for (int k = 0; k < 5; ++k)
        CHECK(r.source_parameters.q[k] == models[i].q[k]);
      if (i < 3) {
        CHECK(r.log_density == e.result.density.log_value &&
              r.quadratic == e.result.quadratic &&
              r.log_determinant == e.result.log_determinant &&
              r.normalization == e.result.normalization);
        CHECK(r.predictions.length == 2 && r.residuals.length == 2);
        for (int k = 0; k < 2; ++k)
          CHECK(r.predictions.data[k] == e.predictions[k] &&
                r.residuals.data[k] == e.residuals[k]);
      } else
        CHECK(r.predictions.length == 0 && r.residuals.length == 0 &&
              r.log_density == 0 && r.quadratic == 0);
    }
    // Analytic independent scalar check only; coarse parity is not science.
    CHECK(std::abs(rows[0].quadratic - .6) < 1e-13);
    CHECK(std::abs(rows[0].log_density +
                   .5 * (.6 + std::log(35.) +
                         2 * std::log(2 * std::acos(-1.)))) < 1e-13);
    f.observed[0] = 999;
    models[0].q[0] = 2;
    CHECK(source.source.observed.data[0] != 999 &&
          rows[0].source_parameters.q[0] == -1);
    CHECK(cosmo_bao_piecewise_destroy(owner) == COSMO_OK);
    owner = nullptr;
    CHECK(string(rows[0].ordered_ids.data[0]) == f.a &&
          string(rows[0].ordered_ids.data[1]) == f.b &&
          rows[0].queries[0].z == .1);
    CHECK(rows[0].predictions.data[0] == expected.slots[0].predictions[0]);
    CHECK(cosmo_bao_piecewise_result_destroy(result) == COSMO_OK);
    result = nullptr;
    f.observed[0] = 299792.458 / 10000 * .1 + 1;
    models[0].q[0] = -1;
    CHECK(cosmo_bao_piecewise_prepare(&f.descriptor, &p, &owner) == COSMO_OK &&
          owner);
    // Prep-only row/matrix/string limits do not reinterpret retained owner.
    auto low = p;
    low.maximum_rows = low.maximum_matrix_elements = low.maximum_string_bytes =
        0;
    CHECK(cosmo_bao_piecewise_evaluate(owner, &batch, &low, &result) ==
              COSMO_OK &&
          result);
    CHECK(cosmo_bao_piecewise_result_view(result, &rows, &count, &status,
                                          &cause, &work) == COSMO_OK &&
          count == 5 && rows[0].status == 0);
    cosmo_bao_piecewise_result_destroy(result);
    result = nullptr;
    for (int k = 0; k < 3; ++k) {
      low = p;
      if (k == 0)
        low.maximum_models = 0;
      if (k == 1) {
        low.maximum_array_elements = 19;
        low.include_predictions = low.include_residuals = 0;
      }
      if (k == 2)
        low.maximum_native_output_bytes = 0;
      CHECK(cosmo_bao_piecewise_evaluate(owner, &batch, &low, &result) ==
                COSMO_OK &&
            result);
      CHECK(cosmo_bao_piecewise_result_view(result, &rows, &count, &status,
                                            &cause, &work) == COSMO_OK &&
            count == 0 && work == 0 &&
            cause == static_cast<uint32_t>(numerics::Status::work_limit));
      cosmo_bao_piecewise_result_destroy(result);
      result = nullptr;
    }
    low = p;
    low.arithmetic = 0;
    CHECK(cosmo_bao_piecewise_evaluate(owner, &batch, &low, &result) ==
              COSMO_OK &&
          result);
    CHECK(cosmo_bao_piecewise_result_view(result, &rows, &count, &status,
                                          &cause, &work) == COSMO_OK &&
          count == 0 && work == 0);
    auto mismatch = native.evaluate_piecewise(points, native_policy(low));
    CHECK(status == static_cast<uint32_t>(mismatch.status) &&
          cause == static_cast<uint32_t>(mismatch.numerical_status));
    cosmo_bao_piecewise_result_destroy(result);
    result = nullptr;
    low = p;
    low.maximum_total_segment_visits = 1;
    CHECK(cosmo_bao_piecewise_evaluate(owner, &batch, &low, &result) ==
              COSMO_OK &&
          result);
    CHECK(cosmo_bao_piecewise_result_view(result, &rows, &count, &status,
                                          &cause, &work) == COSMO_OK &&
          work <= 1);
    auto bounded = native.evaluate_piecewise(points, native_policy(low));
    CHECK(work == bounded.segment_visits && count == bounded.slots.size());
    for (std::size_t i = 0; i < count; ++i)
      CHECK(rows[i].segment_visits == bounded.slots[i].segment_visits &&
            rows[i].status ==
                static_cast<uint32_t>(bounded.slots[i].result.density.status));
    cosmo_bao_piecewise_result_destroy(result);
    result = nullptr;
    auto bad = batch;
    bad.model_count = UINT64_MAX;
    bad.model_byte_length = UINT64_MAX;
    auto before = live;
    CHECK(cosmo_bao_piecewise_evaluate(owner, &bad, &p, &result) ==
              COSMO_INVALID_INPUT &&
          !result && live == before);
    bad = batch;
    --bad.model_byte_length;
    CHECK(cosmo_bao_piecewise_evaluate(owner, &bad, &p, &result) ==
              COSMO_INVALID_INPUT &&
          !result);
    bad = batch;
    bad.abi_version++;
    CHECK(cosmo_bao_piecewise_evaluate(owner, &bad, &p, &result) ==
              COSMO_ABI_MISMATCH &&
          !result);
    low = p;
    low.arithmetic = UINT32_MAX;
    CHECK(cosmo_bao_piecewise_evaluate(owner, &batch, &low, &result) ==
              COSMO_INVALID_INPUT &&
          !result);
    // Complete bad_alloc sweeps: no leaked partial owner/result or unwinding.
    unsigned eval_failures = 0, prep_failures = 0;
    for (int n = 0; n < 512; ++n) {
      before = live;
      fail_after = n;
      auto rc = cosmo_bao_piecewise_evaluate(owner, &batch, &p, &result);
      fail_after = -1;
      if (rc == COSMO_OK) {
        CHECK(result);
        cosmo_bao_piecewise_result_destroy(result);
        result = nullptr;
        CHECK(live == before);
        break;
      }
      CHECK(rc == COSMO_ALLOCATION_FAILURE && !result && live == before);
      ++eval_failures;
    }
    CHECK(eval_failures > 0 && eval_failures < 511);
    // Baseline published archive with this exact fixture admitted legacy
    // preparation at 9745 bytes, rejected 9744. Explicit x86_64 storage check;
    // not a portable floating-point oracle.
    struct LegacyStorageLayout {
      bao::DensityInput source;
      bao::PreparedDensity native;
      cosmo_bao_policy preparation;
      std::vector<cosmo_bao_query> queries;
      std::vector<cosmo_bytes> ids;
      std::string redshift, ruler, h0;
    };
#if defined(__x86_64__) && defined(__GLIBCXX__)
    if (sizeof(LegacyStorageLayout) == 1520 && sizeof(void *) == 8 &&
        sizeof(std::string) == 32 && sizeof(std::vector<double>) == 24 &&
        sizeof(long double) == 16) {
      auto q = f.policy;
      q.maximum_native_bytes = 9744;
      cosmo_bao *old = nullptr;
      CHECK(cosmo_bao_prepare(&f.descriptor, &q, &old) == COSMO_INVALID_INPUT &&
            !old);
      q.maximum_native_bytes = 9745;
      CHECK(cosmo_bao_prepare(&f.descriptor, &q, &old) == COSMO_OK && old);
      CHECK(cosmo_bao_destroy(old) == COSMO_OK);
    }
#endif
    // N=2/P=5 combined metadata, scratch and retained 2Pn arrays. Export
    // suppression must not change admission: native arrays still exist.
    auto admits = [&](uint64_t cap, bool exports) {
      auto q = p;
      q.maximum_native_output_bytes = cap;
      q.include_predictions = q.include_residuals = exports;
      CHECK(cosmo_bao_piecewise_evaluate(owner, &batch, &q, &result) ==
                COSMO_OK &&
            result);
      CHECK(cosmo_bao_piecewise_result_view(result, &rows, &count, &status,
                                            &cause, &work) == COSMO_OK);
      const bool yes = count == 5 && rows[0].status == 0;
      CHECK(yes ||
            (count == 0 && work == 0 &&
             cause == static_cast<uint32_t>(numerics::Status::work_limit)));
      CHECK(cosmo_bao_piecewise_result_destroy(result) == COSMO_OK);
      result = nullptr;
      return yes;
    };
    uint64_t left = 0, right = p.maximum_native_output_bytes;
    while (left < right) {
      auto middle = left + (right - left) / 2;
      if (admits(middle, true))
        right = middle;
      else
        left = middle + 1;
    }
    CHECK(left > 0 && admits(left, true) && !admits(left - 1, true));
    CHECK(admits(left, false) && !admits(left - 1, false));
    cosmo_bao_piecewise_destroy(owner);
    owner = nullptr;
    for (int n = 0; n < 512; ++n) {
      before = live;
      fail_after = n;
      auto rc = cosmo_bao_piecewise_prepare(&f.descriptor, &p, &owner);
      fail_after = -1;
      if (rc == COSMO_OK) {
        CHECK(owner);
        cosmo_bao_piecewise_destroy(owner);
        owner = nullptr;
        CHECK(live == before);
        break;
      }
      CHECK(rc == COSMO_ALLOCATION_FAILURE && !owner && live == before);
      ++prep_failures;
    }
    CHECK(prep_failures > 0 && prep_failures < 511);
    auto descriptor = f.descriptor;
    descriptor.query_count = descriptor.query_byte_length = UINT64_MAX;
    before = live;
    CHECK(cosmo_bao_piecewise_prepare(&descriptor, &p, &owner) ==
              COSMO_INVALID_INPUT &&
          !owner && live == before);
    descriptor = f.descriptor;
    descriptor.abi_version++;
    CHECK(cosmo_bao_piecewise_prepare(&descriptor, &p, &owner) ==
              COSMO_ABI_MISMATCH &&
          !owner);
    descriptor = f.descriptor;
    descriptor.role = UINT32_MAX;
    CHECK(cosmo_bao_piecewise_prepare(&descriptor, &p, &owner) ==
              COSMO_INVALID_INPUT &&
          !owner);
    auto prep_cap = p;
    prep_cap.maximum_matrix_elements = 3;
    CHECK(cosmo_bao_piecewise_prepare(&f.descriptor, &prep_cap, &owner) ==
              COSMO_INVALID_INPUT &&
          !owner);
    // Well-typed unknown observable is retained as an owned scientific
    // preparation failure, not coerced to an existing physical quantity.
    f.queries[0].observable = UINT32_MAX;
    CHECK(cosmo_bao_piecewise_prepare(&f.descriptor, &p, &owner) == COSMO_OK &&
          owner);
    CHECK(cosmo_bao_piecewise_source_view(owner, &source) == COSMO_OK &&
          source.status != COSMO_GAUSSIAN_STATUS_FINITE &&
          source.source.queries[0].observable == UINT32_MAX);
    CHECK(cosmo_bao_piecewise_destroy(owner) == COSMO_OK);
    owner = nullptr;
    std::printf("BAO piecewise ABI independent PASS %d allocationfailures "
                "%u/%u aggregateboundary %llu\n",
                checks, prep_failures, eval_failures,
                static_cast<unsigned long long>(left));
    return 0;
  } catch (const std::exception &e) {
    fail_after = -1;
    std::fprintf(stderr, "BAO piecewise ABI FAIL %s after %d\n", e.what(),
                 checks);
    return 1;
  }
}
