// Independent retained analytic SN transport tests. Generated release-shaped
// controls are not actual survey science. Native/refinement gates are separate.
#include "irred/abi.h"
#include "irred/observations.hpp"
#include "irred/supernova.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <string_view>
static long long fail_after = -1, live = 0;
static unsigned checks = 0;
static void (*volatile release_allocation)(void *) = std::free;
void *operator new(std::size_t n) {
  if (fail_after == 0)
    throw std::bad_alloc();
  if (fail_after > 0)
    --fail_after;
  auto *p = std::malloc(n ? n : 1);
  if (!p)
    throw std::bad_alloc();
  ++live;
  return p;
}
void operator delete(void *p) noexcept {
  if (p) {
    --live;
    release_allocation(p);
  }
}
void operator delete(void *p, std::size_t) noexcept { ::operator delete(p); }
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete[](void *p, std::size_t) noexcept { ::operator delete(p); }
#define CHECK(x)                                                               \
  do {                                                                         \
    ++checks;                                                                  \
    if (!(x))                                                                  \
      throw std::runtime_error("piecewise SN ABI: " #x);                       \
  } while (false)
using namespace irred;
cosmo_bytes bytes(const std::string &s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
std::string_view text(cosmo_bytes b) {
  return {reinterpret_cast<const char *>(b.data),
          static_cast<std::size_t>(b.length)};
}
cosmo_f64_buffer f64(const double *p, std::size_t n) {
  return {sizeof(cosmo_f64_buffer), COSMO_ABI_VERSION, 2, 0, p, n,
          n * sizeof(double)};
}
struct Fixture {
  std::string hash = std::string(64, 'a'), cov = std::string(64, 'b'),
              excluded = "excluded", a = "a", b = "b", event = "same event",
              order =
                  "generated release-shaped boundary control, supplied axes",
              cal = "generated fitted-magnitude control, not actual release",
              dep = "unknown cross-probe dependence";
  cosmo_bytes ids[3] = {bytes(excluded), bytes(a), bytes(b)},
              events[3] = {bytes(event), bytes(event), bytes(event)};
  double value[3] = {100, 2, -3}, matrix[9] = {-1, 7, 8, -7, 4, 1, -8, 1, 9},
         zhd[3] = {.01, 1, 2}, zhel[3] = {.01, .9, 2.1};
  uint8_t missing[3] = {};
  uint64_t quality[3] = {};
  cosmo_observation_policy op{3, 9, 32768};
  cosmo_observation_descriptor d{};
  Fixture() {
    using namespace observations;
    d.struct_size = sizeof(d);
    d.abi_version = COSMO_ABI_VERSION;
    d.profile = static_cast<uint32_t>(Profile::pantheon_plus_released_v1);
    d.role = static_cast<uint32_t>(observations::Role::released_fitted_summary);
    d.unit = static_cast<uint32_t>(observations::Unit::magnitude);
    d.calibration = static_cast<uint32_t>(Calibration::unknown);
    d.uncertainty = static_cast<uint32_t>(Uncertainty::covariance);
    d.uncertainty_unit =
        static_cast<uint32_t>(observations::UncertaintyUnit::magnitude_squared);
    d.table_sha256 = bytes(hash);
    d.uncertainty_sha256 = bytes(cov);
    d.ordering_provenance = bytes(order);
    d.calibration_provenance = bytes(cal);
    d.dependence_provenance = bytes(dep);
    d.measurement_ids = d.uncertainty_axis_ids = {ids, 3, sizeof ids};
    d.event_ids = {events, 3, sizeof events};
    d.values = f64(value, 3);
    d.zhd = d.zcmb = f64(zhd, 3);
    d.zhel = f64(zhel, 3);
    d.uncertainty_matrix = f64(matrix, 9);
    d.missing = d.zhd_missing = d.zcmb_missing =
        d.zhel_missing = {missing, 3, 3};
    d.quality = {quality, 3, sizeof quality};
  }
  observations::Prepared native() {
    using namespace observations;
    Input x;
    x.profile = Profile::pantheon_plus_released_v1;
    x.role = observations::Role::released_fitted_summary;
    x.unit = observations::Unit::magnitude;
    x.calibration = Calibration::unknown;
    x.uncertainty = Uncertainty::covariance;
    x.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
    x.table_sha256 = hash;
    x.uncertainty_sha256 = cov;
    x.ordering_provenance = order;
    x.calibration_provenance = cal;
    x.dependence_provenance = dep;
    x.measurement_ids = x.uncertainty_axis_ids = {excluded, a, b};
    x.event_ids = {event, event, event};
    x.values = {value[0], value[1], value[2]};
    x.zhd = x.zcmb = {zhd[0], zhd[1], zhd[2]};
    x.zhel = {zhel[0], zhel[1], zhel[2]};
    x.missing = x.zhd_missing = x.zcmb_missing = x.zhel_missing = {0, 0, 0};
    x.quality = {0, 0, 0};
    x.uncertainty_matrix.assign(matrix, matrix + 9);
    return prepare(std::move(x), {3, 9, 32768});
  }
};
cosmo_supernova_piecewise_policy policy() {
  return {sizeof(cosmo_supernova_piecewise_policy),
          COSMO_ABI_VERSION,
          1,
          1,
          4,
          3,
          9,
          3,
          64,
          1000000,
          100,
          1e-10};
}
cosmo_supernova_policy legacy() {
  cosmo_supernova_policy p{};
  p.struct_size = sizeof(p);
  p.abi_version = COSMO_ABI_VERSION;
  p.arithmetic = 1;
  p.include_residual_arrays = 1;
  p.maximum_models = 4;
  p.maximum_source_rows = 3;
  p.maximum_matrix_elements = 9;
  p.maximum_array_elements = 64;
  p.maximum_native_output_bytes = 1000000;
  p.maximum_total_evaluations = 100000;
  p.maximum_evaluations_per_integral = 10000;
  p.maximum_depth = 24;
  p.absolute_tolerance = p.relative_tolerance = 1e-12;
  p.maximum_forward_sensitivity = 1e-10;
  return p;
}
cosmo_prepared *observations_owner(Fixture &f) {
  cosmo_prepared *r = nullptr;
  uint32_t status = 999;
  CHECK(cosmo_prepare_observations(&f.d, &f.op, &r, &status) == COSMO_OK && r &&
        status == 0);
  return r;
}
struct View {
  const cosmo_supernova_piecewise_slot *rows = nullptr;
  uint64_t count = 0, visits = 0;
  uint32_t status = 999;
};
View view(cosmo_supernova_piecewise_result *r) {
  View v;
  CHECK(cosmo_supernova_piecewise_result_view(r, &v.rows, &v.count, &v.status,
                                              &v.visits) == COSMO_OK);
  return v;
}
void absent(const cosmo_supernova_piecewise_slot &r) {
  CHECK(r.has_profile_payload == 0);
  CHECK(r.shape_magnitudes.length == 0 && r.base_residuals.length == 0 &&
        r.profiled_residuals.length == 0);
  CHECK(r.relative_profile_score == 0 && r.quadratic == 0 &&
        r.offset_coefficient == 0);
}
int main() {
  try {
    static_assert(sizeof(cosmo_supernova_piecewise_model) == 40 &&
                  sizeof(cosmo_supernova_piecewise_batch) == 32 &&
                  sizeof(cosmo_supernova_piecewise_policy) == 80 &&
                  sizeof(cosmo_supernova_piecewise_slot) == 488);
    static_assert(sizeof(cosmo_supernova_model) == 24 &&
                  sizeof(cosmo_supernova_model_v2) == 40);
    CHECK(offsetof(cosmo_supernova_piecewise_slot, source_parameters) == 16);
    Fixture f;
    auto *obs = observations_owner(f);
    auto p = policy();
    cosmo_supernova *owner = nullptr, *old = nullptr;
    CHECK(cosmo_supernova_piecewise_prepare(obs, &p, &owner) == COSMO_OK &&
          owner);
    auto lp = legacy();
    CHECK(cosmo_supernova_prepare(obs, &lp, &old) == COSMO_OK && old);
    cosmo_supernova_view sv{}, ov{};
    CHECK(cosmo_supernova_source_view(owner, &sv) == COSMO_OK &&
          sv.status == 0);
    CHECK(cosmo_supernova_source_view(old, &ov) == COSMO_OK &&
          ov.status == sv.status);
    CHECK(sv.source_row_count == 3 && sv.ordered_ids.length == 2 &&
          sv.selected_source_indices.data[0] == 1 &&
          sv.selected_source_indices.data[1] == 2);
    CHECK(sv.matrix_validation_assessed == 1 &&
          sv.matrix_validation_scope == ov.matrix_validation_scope);
    CHECK(text(sv.ordering_provenance) == f.order &&
          text(sv.calibration_provenance) == f.cal &&
          text(sv.dependence_provenance) == f.dep);
    CHECK(sv.selected_z_expansion.data[0] == 1 &&
          sv.selected_z_observer.data[1] == 2.1);
    // Old prepare full-source n^2 guard and scientific failure view are
    // unchanged.
    for (bool analytic : {false, true}) {
      cosmo_supernova *tmp = nullptr;
      if (analytic) {
        auto bad = p;
        bad.maximum_matrix_elements = 8;
        CHECK(cosmo_supernova_piecewise_prepare(obs, &bad, &tmp) ==
                  COSMO_INVALID_INPUT &&
              !tmp);
      } else {
        auto bad = lp;
        bad.maximum_matrix_elements = 8;
        CHECK(cosmo_supernova_prepare(obs, &bad, &tmp) == COSMO_INVALID_INPUT &&
              !tmp);
      }
    }
    auto tight = p;
    tight.maximum_forward_sensitivity = 1e-30;
    cosmo_supernova *failed = nullptr;
    CHECK(cosmo_supernova_piecewise_prepare(obs, &tight, &failed) == COSMO_OK &&
          failed);
    cosmo_supernova_view fv{};
    CHECK(cosmo_supernova_source_view(failed, &fv) == COSMO_OK &&
          fv.status != 0);
    CHECK(fv.preparation_numerical_status ==
              COSMO_NUMERICAL_STATUS_CONDITIONING_BUDGET_EXCEEDED &&
          fv.matrix_validation_assessed == 0 && fv.arithmetic_id.length == 0);
    CHECK(fv.source_row_count == 3 && text(fv.table_identity) == f.hash);
    cosmo_supernova_destroy(failed);
    // Actual prepare/evaluate allocations are all exercised, without a fixed
    // count.
    unsigned prep_failures = 0, eval_failures = 0;
    bool reached = false;
    for (int k = 0; k < 300; ++k) {
      auto base = live;
      cosmo_supernova *tmp = nullptr;
      fail_after = k;
      auto st = cosmo_supernova_piecewise_prepare(obs, &p, &tmp);
      fail_after = -1;
      if (st == COSMO_OK) {
        CHECK(tmp);
        cosmo_supernova_destroy(tmp);
        reached = true;
      } else {
        CHECK(st == COSMO_ALLOCATION_FAILURE && !tmp);
        ++prep_failures;
      }
      CHECK(live == base);
      if (reached)
        break;
    }
    CHECK(reached && prep_failures > 0);
    supernova::Policy np;
    np.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
    np.maximum_forward_sensitivity = 1e-10;
    np.maximum_matrix_elements = 9;
    np.background.maximum_queries = 3;
    auto direct = supernova::prepare(f.native(), np);
    CHECK(direct.status() == supernova::Status::ok);
    std::array<cosmo_supernova_piecewise_model, 3> models{
        {{{-1, -1, -1, -1, -1}}, {{0, 0, 0, 0, 0}}, {{3, 0, 0, 0, 0}}}};
    cosmo_supernova_piecewise_batch batch{sizeof(batch), COSMO_ABI_VERSION,
                                          models.data(), models.size(),
                                          sizeof(models)};
    std::array<supernova::PiecewiseModelPoint, 3> points{
        supernova::PiecewiseModelPoint({-1, -1, -1, -1, -1}),
        supernova::PiecewiseModelPoint({0, 0, 0, 0, 0}),
        supernova::PiecewiseModelPoint({3, 0, 0, 0, 0})};
    supernova::PiecewiseEvaluationPolicy ep;
    ep.arithmetic = np.arithmetic;
    ep.maximum_forward_sensitivity = 1e-10;
    ep.background.maximum_queries = 3;
    ep.background.maximum_segment_visits = 100;
    auto nb = direct.evaluate_piecewise_batch(points, ep);
    cosmo_supernova_piecewise_result *r = nullptr;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &p, &r) ==
              COSMO_OK &&
          r);
    auto v = view(r);
    CHECK(v.count == 3 && v.status == 0 && v.visits == nb.segment_visits);
    for (std::size_t i = 0; i < 3; ++i) {
      const auto &row = v.rows[i];
      const auto &s = nb.slots[i];
      CHECK(row.model_index == i && row.struct_size == sizeof(row) &&
            row.abi_version == COSMO_ABI_VERSION);
      CHECK(row.status == static_cast<uint32_t>(s.status) &&
            row.background_status ==
                static_cast<uint32_t>(s.background_status) &&
            row.numerical_status == static_cast<uint32_t>(s.numerical_status));
      for (int k = 0; k < 5; ++k)
        CHECK(std::bit_cast<uint64_t>(row.source_parameters.q[k]) ==
              std::bit_cast<uint64_t>(models[i].q[k]));
      CHECK(row.segment_visits == s.segment_visits);
      CHECK(row.ordered_ids.length == 2 && row.selected_count == 2 &&
            row.selected_source_indices[0] == 1 &&
            row.selected_source_indices[1] == 2);
      CHECK(row.expansion_z.data[0] == 1 && row.observer_z.data[1] == 2.1);
      if (s.status != supernova::Status::ok) {
        absent(row);
        continue;
      }
      CHECK(row.has_profile_payload == 1 &&
            row.relative_profile_score == s.relative_profile_score &&
            row.quadratic == s.quadratic &&
            row.offset_coefficient == s.offset_coefficient);
      CHECK(row.backward_residual == s.solve_diagnostics.backward_residual &&
            row.estimated_forward_sensitivity ==
                s.solve_diagnostics.estimated_forward_sensitivity);
      for (std::size_t j = 0; j < 2; ++j)
        CHECK(row.shape_magnitudes.data[j] == s.shape_magnitudes[j] &&
              row.base_residuals.data[j] == s.base_residuals[j] &&
              row.profiled_residuals.data[j] ==
                  s.solve_diagnostics.adjusted_residuals[j]);
    }
    long double r0 = 2 - 5 * std::log10(1.9L), r1 = -3 - 5 * std::log10(6.2L);
    CHECK(std::abs(v.rows[0].quadratic -
                   static_cast<double>((r0 - r1) * (r0 - r1) / 11)) < 1e-10);
    cosmo_supernova_piecewise_result_destroy(r);
    reached = false;
    for (int k = 0; k < 300; ++k) {
      auto base = live;
      r = nullptr;
      fail_after = k;
      auto st = cosmo_supernova_piecewise_evaluate(owner, &batch, &p, &r);
      fail_after = -1;
      if (st == COSMO_OK) {
        CHECK(r);
        cosmo_supernova_piecewise_result_destroy(r);
        reached = true;
      } else {
        CHECK(st == COSMO_ALLOCATION_FAILURE && !r);
        ++eval_failures;
      }
      CHECK(live == base);
      if (reached)
        break;
    }
    CHECK(reached && eval_failures > 0);
    // Export suppression still accounts four actual native arrays per P*n.
    auto silent = p;
    silent.include_residual_arrays = 0;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &silent, &r) ==
              COSMO_OK &&
          r);
    v = view(r);
    CHECK(v.rows[0].has_profile_payload == 1 &&
          v.rows[0].shape_magnitudes.length == 0);
    cosmo_supernova_piecewise_result_destroy(r);
    silent.maximum_array_elements = 23;
    r = nullptr;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &silent, &r) ==
              COSMO_OK &&
          r);
    v = view(r);
    CHECK(v.status == COSMO_SUPERNOVA_STATUS_WORK_LIMIT && v.count == 0 &&
          v.visits == 0);
    cosmo_supernova_piecewise_result_destroy(r);
    // Preparation dimensions do not become selected-row evaluation caps.
    auto phase = p;
    phase.maximum_source_rows = 0;
    phase.maximum_matrix_elements = 0;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &phase, &r) ==
              COSMO_OK && r);
    v = view(r);
    CHECK(v.rows[0].has_profile_payload == 1 &&
          v.rows[0].relative_profile_score == nb.slots[0].relative_profile_score);
    cosmo_supernova_piecewise_result_destroy(r);
    auto mismatch = p;
    mismatch.arithmetic = 0;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &mismatch, &r) ==
              COSMO_OK &&
          r);
    v = view(r);
    CHECK(v.status == COSMO_SUPERNOVA_STATUS_INCOMPATIBLE_METADATA &&
          v.count == 0);
    cosmo_supernova_piecewise_result_destroy(r);
    tight = p;
    tight.maximum_forward_sensitivity = 1e-30;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &tight, &r) ==
              COSMO_OK &&
          r);
    v = view(r);
    CHECK(v.rows[0].numerical_status ==
              COSMO_NUMERICAL_STATUS_CONDITIONING_BUDGET_EXCEEDED &&
          v.rows[0].segment_visits == 9);
    absent(v.rows[0]);
    cosmo_supernova_piecewise_result_destroy(r);
    auto limited = p;
    limited.maximum_total_segment_visits = 9;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &limited, &r) ==
              COSMO_OK &&
          r);
    v = view(r);
    CHECK(v.visits == 9 && v.rows[0].has_profile_payload == 1 &&
          v.rows[1].background_status == COSMO_BACKGROUND_STATUS_WORK_LIMIT &&
          v.rows[1].segment_visits == 0);
    absent(v.rows[1]);
    cosmo_supernova_piecewise_result_destroy(r);
    // Malformed storage/caps must stop before dereference/allocation.
    for (int mode = 0; mode < 5; ++mode) {
      auto b = batch;
      auto bad = p;
      if (mode == 0)
        b.abi_version++;
      if (mode == 1)
        bad.struct_size--;
      if (mode == 2)
        b.model_bytes--;
      if (mode == 3)
        b.models = reinterpret_cast<const cosmo_supernova_piecewise_model *>(
            reinterpret_cast<const char *>(models.data()) + 1);
      if (mode == 4) {
        b.model_count = UINT64_MAX;
        b.model_bytes = UINT64_MAX;
      }
      r = reinterpret_cast<cosmo_supernova_piecewise_result *>(uintptr_t(1));
      auto st = cosmo_supernova_piecewise_evaluate(owner, &b, &bad, &r);
      CHECK(st == (mode == 0 ? COSMO_ABI_MISMATCH : COSMO_INVALID_INPUT) && !r);
    }
    // New valid resource budgets return an owned empty diagnostic. Zero-product
    // admitted batches are OK; malformed policy/shape remains null above.
    for (int mode = 0; mode < 3; ++mode) {
      auto bad = p;
      if (mode == 0)
        bad.maximum_models = 0;
      if (mode == 1)
        bad.maximum_native_output_bytes = 0;
      if (mode == 2)
        bad.maximum_native_output_bytes = 1;
      r = nullptr;
      CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &bad, &r) ==
                COSMO_OK &&
            r);
      v = view(r);
      CHECK(v.status == COSMO_SUPERNOVA_STATUS_WORK_LIMIT && v.count == 0 &&
            v.visits == 0);
      cosmo_supernova_piecewise_result_destroy(r);
    }
    auto empty = batch;
    empty.models = nullptr;
    empty.model_count = 0;
    empty.model_bytes = 0;
    auto zero_models = p;
    zero_models.maximum_models = 0;
    r = nullptr;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &empty, &zero_models, &r) ==
              COSMO_OK &&
          r);
    v = view(r);
    CHECK(v.status == COSMO_SUPERNOVA_STATUS_OK && v.count == 0 &&
          v.visits == 0);
    cosmo_supernova_piecewise_result_destroy(r);
    // Selected-query cap is a native per-row cause, distinct from full source
    // preparation bounds and call-level output caps.
    auto query_limited = p;
    query_limited.maximum_queries = 1;
    CHECK(cosmo_supernova_piecewise_evaluate(owner, &batch, &query_limited,
                                             &r) == COSMO_OK &&
          r);
    v = view(r);
    CHECK(v.count == 3 &&
          v.rows[0].background_status == COSMO_BACKGROUND_STATUS_WORK_LIMIT &&
          v.rows[0].numerical_status == COSMO_NUMERICAL_STATUS_WORK_LIMIT &&
          v.visits == 0);
    absent(v.rows[0]);
    cosmo_supernova_piecewise_result_destroy(r);
    // Old and new entry points share one private owner: old evaluations remain
    // valid.
    cosmo_supernova_model oldpoint{1, 0, 0, -1};
    cosmo_supernova_batch ob{sizeof(ob), COSMO_ABI_VERSION, &oldpoint, 1,
                             sizeof(oldpoint)};
    cosmo_supernova_result *orow = nullptr;
    CHECK(cosmo_supernova_evaluate(owner, &ob, &lp, &orow) == COSMO_OK && orow);
    const cosmo_supernova_slot *oldrows = nullptr;
    uint64_t count = 0;
    uint32_t status = 9;
    CHECK(cosmo_supernova_result_view(orow, &oldrows, &count, &status) ==
              COSMO_OK &&
          count == 1 && status == 0);
    CHECK(std::abs(oldrows[0].relative_profile_score -
                   nb.slots[0].relative_profile_score) < 1e-10);
    cosmo_supernova_result_destroy(orow);
    cosmo_supernova_destroy(old);
    cosmo_supernova_destroy(owner);
    cosmo_observation_destroy(obs);
    // Exercise short SSO and long source IDs after both source owners are
    // destroyed.
    for (bool long_id : {false, true}) {
      Fixture lf;
      if (long_id) {
        lf.a = std::string(1024, 'x');
        lf.ids[1] = bytes(lf.a);
      }
      auto *o = observations_owner(lf);
      cosmo_supernova *h = nullptr;
      CHECK(cosmo_supernova_piecewise_prepare(o, &p, &h) == COSMO_OK && h);
      CHECK(cosmo_supernova_piecewise_evaluate(h, &batch, &p, &r) == COSMO_OK &&
            r);
      const auto expected = lf.a;
      cosmo_supernova_destroy(h);
      cosmo_observation_destroy(o);
      std::vector<std::string> churn(100, std::string(1024, '!'));
      v = view(r);
      CHECK(text(v.rows[0].ordered_ids.data[0]) == expected &&
            text(v.rows[0].ordered_ids.data[1]) == "b");
      CHECK(v.rows[0].selected_source_indices[0] == 1 &&
            v.rows[0].expansion_z.data[0] == 1 &&
            v.rows[0].observer_z.data[1] == 2.1);
      CHECK(v.rows[0].shape_magnitudes.length == 2 &&
            std::isfinite(v.rows[0].shape_magnitudes.data[0]));
      cosmo_supernova_piecewise_result_destroy(r);
    }
    CHECK(cosmo_supernova_piecewise_result_destroy(nullptr) == COSMO_OK);
    std::printf("piecewise SN ABI checks %u PASS allocation failures prepare "
                "%u evaluate %u\n",
                checks, prep_failures, eval_failures);
    return 0;
  } catch (const std::exception &e) {
    fail_after = -1;
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
