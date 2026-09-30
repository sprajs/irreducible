// Original bridge/status/ownership controls authored by statistics core owner;
// current ABI2 cause/selection controls added by the independent verifier. Direct
// parity is NOT independent scientific evidence. Native scientific reference
// qualification lives in test_statistics_hostile. Structural descriptor fixture
// follows test_observation_abi_hostile; no Rust/shared equation duplication.
#include "irred/abi.h"
#include "irred/observations.hpp"
#include "irred/statistics.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
static std::int64_t fail_after = -1, live = 0;
static int checks = 0;
struct alignas(std::max_align_t) AllocationHeader { size_t bytes; bool tracked; };
static bool track_payload=false;
static size_t payload=0,payload_peak=0;
void *operator new(std::size_t n) {
  if (fail_after == 0)
    throw std::bad_alloc();
  if (fail_after > 0)
    --fail_after;
  void *p = std::malloc(sizeof(AllocationHeader)+(n ? n : 1));
  if (!p)
    throw std::bad_alloc();
  ++live;
  auto h=new(p) AllocationHeader{n,track_payload};
  if(track_payload) { payload+=n;payload_peak=std::max(payload_peak,payload); }
  return h+1;
}
void operator delete(void *p) noexcept {
  if (p) {
    --live;
    auto h=static_cast<AllocationHeader*>(p)-1;
    if(h->tracked) payload-=h->bytes;
    std::free(h);
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
      throw std::runtime_error("Gaussian ABI: " #__VA_ARGS__);                 \
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
  double values[2] = {2, -3}, matrix[4] = {4, 1, 1, 9};
  uint8_t missing[2] = {0, 0};
  uint64_t quality[2] = {0, 0};
  cosmo_observation_policy op{2, 4, 4096};
  cosmo_observation_descriptor d{};
  Fixture() {
    d.struct_size = sizeof(d);
    d.abi_version = COSMO_ABI_VERSION;
    d.profile = static_cast<uint32_t>(Profile::gaussian_fixture_v1);
    d.role = static_cast<uint32_t>(Role::synthetic_control);
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
    d.zhd = d.zcmb = d.zhel = doubles(nullptr, 0);
    d.uncertainty_matrix = doubles(matrix, 4);
    d.missing = {missing, 2, 2};
    d.quality = {quality, 2, sizeof quality};
  }
};
int main() {
  try {
    Fixture f;
    cosmo_prepared *obs = nullptr;
    uint32_t semantic = 99, numerical = 99;
    CHECK(cosmo_prepare_observations(&f.d, &f.op, &obs, &semantic) ==
              COSMO_OK &&
          obs);
    cosmo_gaussian_policy policy{sizeof(cosmo_gaussian_policy),
                                 COSMO_ABI_VERSION,
                                 0,
                                 0,
                                 4,
                                 100,
                                 4096,
                                 1e-10,
                                 1000000};
    cosmo_gaussian *g = nullptr;
    CHECK(cosmo_gaussian_prepare(obs, 0, &policy, &g, &semantic, &numerical) == COSMO_OK &&
          g && semantic == COSMO_GAUSSIAN_STATUS_FINITE &&
          numerical == (uint32_t)irred::numerics::Status::ok);
    // Finite typed matrices are structural observations even when Gaussian
    // preparation cannot admit SPD. Preserve transport/semantic/cause layers.
    for (double leading : {0., -1.}) {
      Fixture invalid;
      invalid.matrix[0] = leading;
      invalid.matrix[1] = invalid.matrix[2] = 0;
      cosmo_prepared *source = nullptr;
      CHECK(cosmo_prepare_observations(&invalid.d, &invalid.op, &source,
                                       &semantic) == COSMO_OK && source);
      cosmo_gaussian *failed = nullptr;
      CHECK(cosmo_gaussian_prepare(source, 0, &policy, &failed, &semantic,
                                   &numerical) == COSMO_OK && !failed &&
            semantic == COSMO_GAUSSIAN_STATUS_NUMERICAL_FAILURE &&
            numerical == (uint32_t)irred::numerics::Status::not_positive_definite);
      CHECK(cosmo_observation_destroy(source) == COSMO_OK);
    }
    auto tight_prepare = policy;
    tight_prepare.maximum_forward_sensitivity = 1e-20;
    Fixture precision_source;
    precision_source.d.uncertainty = (uint32_t)Uncertainty::precision;
    precision_source.d.uncertainty_unit =
        (uint32_t)UncertaintyUnit::inverse_magnitude_squared;
    cosmo_prepared *precision_observations = nullptr;
    CHECK(cosmo_prepare_observations(&precision_source.d, &precision_source.op,
                                     &precision_observations, &semantic) == COSMO_OK &&
          precision_observations);
    cosmo_gaussian *tight_failed = nullptr;
    CHECK(cosmo_gaussian_prepare(precision_observations, 0, &tight_prepare, &tight_failed,
                                 &semantic, &numerical) == COSMO_OK &&
          !tight_failed && semantic == COSMO_GAUSSIAN_STATUS_NUMERICAL_FAILURE &&
          numerical == (uint32_t)irred::numerics::Status::conditioning_budget_exceeded);
    CHECK(cosmo_observation_destroy(precision_observations) == COSMO_OK);
    for (bool precision : {false,true}) {
      Fixture selected;
      uint8_t selected_mask[2]{1,0};
      selected.d.source_selection = {selected_mask,2,2};
      selected.d.uncertainty = (uint32_t)(precision ? Uncertainty::precision
                                                   : Uncertainty::covariance);
      selected.d.uncertainty_unit = (uint32_t)(precision
          ? UncertaintyUnit::inverse_magnitude_squared
          : UncertaintyUnit::magnitude_squared);
      if (!precision) selected.matrix[3] = -9;
      cosmo_prepared *source = nullptr;
      CHECK(cosmo_prepare_observations(&selected.d, &selected.op, &source,
                                       &semantic) == COSMO_OK && source);
      cosmo_gaussian *marginal = nullptr;
      CHECK(cosmo_gaussian_prepare(source,0,&policy,&marginal,&semantic,
                                   &numerical) == COSMO_OK && marginal &&
            semantic == COSMO_GAUSSIAN_STATUS_FINITE &&
            numerical == (uint32_t)irred::numerics::Status::ok);
      cosmo_gaussian_view source_view{};
      CHECK(cosmo_gaussian_source_view(marginal,&source_view) == COSMO_OK &&
            source_view.ordered_ids.length == 1 &&
            source_view.ordered_ids.data[0].data[0] == 'a' &&
            source_view.matrix_validation_scope == (uint32_t)(precision
              ? irred::statistics::MatrixValidationScope::full_precision_then_marginal
              : irred::statistics::MatrixValidationScope::selected_covariance_only));
      double residual = 2;
      cosmo_gaussian_batch one{sizeof(one),COSMO_ABI_VERSION,
          COSMO_GAUSSIAN_MODE_NORMALIZED_DENSITY,0,1,doubles(&residual,1),
          {selected.ids,1,sizeof(cosmo_bytes)},doubles(nullptr,0)};
      cosmo_gaussian_result *result = nullptr;
      CHECK(cosmo_gaussian_evaluate(marginal,&one,&policy,&result) == COSMO_OK && result);
      const cosmo_gaussian_row *row = nullptr;
      uint64_t count = 0;
      const double variance = precision ? 9./35 : 4.;
      CHECK(cosmo_gaussian_result_view(result,&row,&count) == COSMO_OK && count == 1 &&
            row[0].status == COSMO_GAUSSIAN_STATUS_FINITE &&
            std::abs(row[0].quadratic-4/variance) < 4e-11 &&
            std::abs(row[0].log_determinant-std::log(variance)) < 4e-11 &&
            std::abs(row[0].normalization-std::log(2*std::acos(-1.))) < 4e-11);
      CHECK(cosmo_gaussian_result_destroy(result) == COSMO_OK);
      CHECK(cosmo_gaussian_destroy(marginal) == COSMO_OK);
      CHECK(cosmo_observation_destroy(source) == COSMO_OK);
    }
    bool prepared_success = false;
    int prepare_failures = 0;
    for (int point = 0; point < 300; ++point) {
      auto before = live;
      cosmo_gaussian *temporary = nullptr;
      fail_after = point;
      auto status =
          cosmo_gaussian_prepare(obs, 0, &policy, &temporary, &semantic, &numerical);
      fail_after = -1;
      if (status == COSMO_OK) {
        CHECK(temporary);
        CHECK(cosmo_gaussian_destroy(temporary) == COSMO_OK);
        CHECK(live == before);
        prepared_success = true;
        break;
      }
      CHECK(status == COSMO_ALLOCATION_FAILURE && !temporary && live == before);
      ++prepare_failures;
    }
    CHECK(prepared_success && prepare_failures > 0);
    auto invalid_policy = policy;
    invalid_policy.abi_version++;
    cosmo_gaussian *temporary = nullptr;
    CHECK(cosmo_gaussian_prepare(obs, 0, &invalid_policy, &temporary,
                                 &semantic, &numerical) == COSMO_ABI_MISMATCH &&
          !temporary);
    CHECK(cosmo_gaussian_prepare(nullptr, 0, &policy, &temporary, &semantic, &numerical) ==
              COSMO_INVALID_INPUT &&
          !temporary);
    CHECK(cosmo_gaussian_prepare(obs, 99, &policy, &temporary, &semantic, &numerical) ==
              COSMO_OK &&
          !temporary);
    CHECK(cosmo_observation_destroy(obs) == COSMO_OK);
    obs = nullptr;
    f.matrix[0] = 123;
    cosmo_gaussian_view view{};
    CHECK(cosmo_gaussian_source_view(g, &view) == COSMO_OK &&
          view.ordered_ids.length == 2 &&
          view.calibration_provenance.length == 0 &&
          view.dependence_provenance.length == 0 && view.prior_count == 0);
    CHECK(view.matrix_validation_scope == 1);
    CHECK(view.selection_count == 1);
    cosmo_gaussian_selection history{};
    CHECK(cosmo_gaussian_selection_view(g, 0, &history) == COSMO_OK &&
          history.kept_row_ids.length == 2 &&
          history.complement_row_ids.length == 0);
    CHECK(cosmo_gaussian_selection_view(g, 1, &history) ==
              COSMO_INVALID_INPUT &&
          history.kept_row_ids.data == nullptr);
    CHECK(std::string(reinterpret_cast<const char *>(view.table_identity.data),
                      view.table_identity.length) == f.hash);
    std::array<double, 6> residual{2, -3, 0, 0, NAN, 1};
    cosmo_gaussian_batch batch{sizeof(cosmo_gaussian_batch),
                               COSMO_ABI_VERSION,
                               COSMO_GAUSSIAN_MODE_NORMALIZED_DENSITY,
                               0,
                               3,
                               doubles(residual.data(), 6),
                               {f.ids, 2, sizeof f.ids},
                               doubles(nullptr, 0)};
    cosmo_gaussian_result *out = nullptr;
    CHECK(cosmo_gaussian_evaluate(g, &batch, &policy, &out) == COSMO_OK && out);
    const cosmo_gaussian_row *rows = nullptr;
    uint64_t n = 0;
    CHECK(cosmo_gaussian_result_view(out, &rows, &n) == COSMO_OK && n == 3);
    CHECK(rows[0].status == COSMO_GAUSSIAN_STATUS_FINITE &&
          std::abs(rows[0].quadratic - 2.4) < 4e-11 &&
          std::abs(rows[0].log_determinant - std::log(35.)) < 4e-11);
    CHECK(rows[1].status == COSMO_GAUSSIAN_STATUS_FINITE &&
          rows[1].quadratic == 0);
    CHECK(rows[2].status == COSMO_GAUSSIAN_STATUS_NUMERICAL_FAILURE);
    CHECK(cosmo_gaussian_result_destroy(out) == COSMO_OK);
    out = nullptr;
    std::array<double, 2> response{1, 1};
    batch.mode = COSMO_GAUSSIAN_MODE_PROFILE_OFFSET_SCORE;
    batch.row_count = 1;
    batch.residuals = doubles(residual.data(), 2);
    batch.response = doubles(response.data(), 2);
    CHECK(cosmo_gaussian_evaluate(g, &batch, &policy, &out) == COSMO_OK && out);
    CHECK(cosmo_gaussian_result_view(out, &rows, &n) == COSMO_OK && n == 1);
    CHECK(rows[0].status == COSMO_GAUSSIAN_STATUS_FINITE &&
          std::abs(rows[0].coefficient - 7. / 11) < 1e-10 &&
          std::abs(rows[0].quadratic - 25. / 11) < 1e-10);
    CHECK(rows[0].numerical_status == COSMO_NUMERICAL_STATUS_OK &&
          rows[0].backward_residual >= 0 &&
          rows[0].estimated_forward_sensitivity > 0);
    CHECK(cosmo_gaussian_result_destroy(out) == COSMO_OK);
    out = nullptr;
    auto tight_policy = policy;
    tight_policy.maximum_forward_sensitivity = 1e-30;
    CHECK(cosmo_gaussian_evaluate(g, &batch, &tight_policy, &out) == COSMO_OK &&
          out);
    CHECK(cosmo_gaussian_result_view(out, &rows, &n) == COSMO_OK && n == 1 &&
          rows[0].status == COSMO_GAUSSIAN_STATUS_NUMERICAL_FAILURE &&
          rows[0].numerical_status ==
              COSMO_NUMERICAL_STATUS_CONDITIONING_BUDGET_EXCEEDED);
    CHECK(cosmo_gaussian_result_destroy(out) == COSMO_OK);
    out = nullptr;
    cosmo_bytes reversed[2] = {f.ids[1], f.ids[0]};
    auto wrong = batch;
    wrong.ordered_ids = {reversed, 2, sizeof reversed};
    CHECK(cosmo_gaussian_evaluate(g, &wrong, &policy, &out) == COSMO_OK && out);
    CHECK(cosmo_gaussian_result_view(out, &rows, &n) == COSMO_OK &&
          rows[0].status == COSMO_GAUSSIAN_STATUS_INCOMPATIBLE_METADATA);
    CHECK(cosmo_gaussian_result_destroy(out) == COSMO_OK);
    out = nullptr;
    cosmo_gaussian_prior prior{sizeof(cosmo_gaussian_prior),
                               COSMO_ABI_VERSION,
                               1,
                               0,
                               1. / 3,
                               2,
                               bytes(f.latent),
                               doubles(response.data(), 2),
                               {f.ids, 2, sizeof f.ids}};
    cosmo_gaussian *proper = nullptr;
    CHECK(cosmo_gaussian_proper_offset(g, &prior, &policy, &proper,
                                       &semantic, &numerical) == COSMO_OK &&
          proper && semantic == COSMO_GAUSSIAN_STATUS_FINITE);
    cosmo_gaussian_prior pv{};
    CHECK(cosmo_gaussian_prior_view(proper, 0, &pv) == COSMO_OK &&
          pv.mean == prior.mean && pv.variance == 2 &&
          pv.independence_declared == 1 && pv.response.data[0] == 1 &&
          pv.ordered_ids.length == 2);
    CHECK(cosmo_gaussian_source_view(proper, &view) == COSMO_OK &&
          view.prior_count == 1 && view.mean_shift.length == 2 &&
          view.mean_shift.data[0] == 1. / 3 &&
          view.dependence_provenance.length == 0);
    CHECK(cosmo_gaussian_evaluate(proper, &batch, &policy, &out) == COSMO_OK);
    CHECK(cosmo_gaussian_result_view(out, &rows, &n) == COSMO_OK &&
          rows[0].status == COSMO_GAUSSIAN_STATUS_INCOMPATIBLE_METADATA);
    CHECK(cosmo_gaussian_result_destroy(out) == COSMO_OK);
    out = nullptr;
    cosmo_gaussian *badg = nullptr;
    CHECK(cosmo_gaussian_proper_offset(proper, &prior, &policy, &badg,
                                       &semantic, &numerical) == COSMO_OK &&
          !badg && semantic != COSMO_GAUSSIAN_STATUS_FINITE);
    auto wrongprior = prior;
    wrongprior.ordered_ids = {reversed, 2, sizeof reversed};
    CHECK(cosmo_gaussian_proper_offset(g, &wrongprior, &policy, &badg,
                                       &semantic, &numerical) == COSMO_OK &&
          !badg && semantic == COSMO_GAUSSIAN_STATUS_INCOMPATIBLE_METADATA);
    auto normalized = batch;
    normalized.mode = COSMO_GAUSSIAN_MODE_NORMALIZED_DENSITY;
    normalized.response = doubles(nullptr, 0);
    CHECK(cosmo_gaussian_evaluate(proper, &normalized, &policy, &out) ==
          COSMO_OK);
    CHECK(cosmo_gaussian_result_view(out, &rows, &n) == COSMO_OK &&
          rows[0].status == COSMO_GAUSSIAN_STATUS_FINITE &&
          std::abs(rows[0].quadratic - 1175. / 513) < 4e-11 &&
          std::abs(rows[0].log_determinant - std::log(57.)) < 4e-11);
    CHECK(cosmo_gaussian_result_destroy(out) == COSMO_OK);
    out = nullptr;
    // Independently measure prepared native+wire owner and result/scratch
    // overlap. Threshold discovery tests admission; observed allocations test
    // that the declared payload envelope actually covers simultaneous storage.
    cosmo_prepared *measurement_source=nullptr;
    CHECK(cosmo_prepare_observations(&f.d,&f.op,&measurement_source,&semantic)==COSMO_OK && measurement_source);
    for (bool profile : {false,true}) {
      CHECK(payload==0);
      track_payload=true;
      cosmo_gaussian *measured=nullptr;
      CHECK(cosmo_gaussian_prepare(measurement_source,0,&policy,&measured,&semantic,&numerical)==COSMO_OK && measured);
      track_payload=false;
      const auto retained_payload=payload;
      auto probe=profile ? batch : normalized;
      uint64_t lo=0,hi=policy.maximum_native_bytes;
      while(lo<hi) {
        auto limit=policy;limit.maximum_native_bytes=lo+(hi-lo)/2;
        cosmo_gaussian_result *candidate=nullptr;
        auto code=cosmo_gaussian_evaluate(measured,&probe,&limit,&candidate);
        if(code==COSMO_OK) {cosmo_gaussian_result_destroy(candidate);hi=limit.maximum_native_bytes;}
        else {CHECK(code==COSMO_INVALID_INPUT && !candidate);lo=limit.maximum_native_bytes+1;}
      }
      auto exact=policy;exact.maximum_native_bytes=lo;
      payload_peak=payload;track_payload=true;
      CHECK(cosmo_gaussian_evaluate(measured,&probe,&exact,&out)==COSMO_OK && out);
      track_payload=false;
      CHECK(payload_peak<=lo);
      CHECK(cosmo_gaussian_result_destroy(out)==COSMO_OK);out=nullptr;
      CHECK(payload==retained_payload);
      --exact.maximum_native_bytes;
      CHECK(cosmo_gaussian_evaluate(measured,&probe,&exact,&out)==COSMO_INVALID_INPUT && !out);
      CHECK(payload==retained_payload);
      std::printf("Gaussian wire evaluation profile %d observed %zu admitted_bound %llu\n",profile,payload_peak,(unsigned long long)lo);
      CHECK(cosmo_gaussian_destroy(measured)==COSMO_OK && payload==0);
    }
    CHECK(cosmo_observation_destroy(measurement_source)==COSMO_OK);
    auto before_views = live;
    fail_after = 0;
    auto source_status = cosmo_gaussian_source_view(proper, &view);
    auto prior_status = cosmo_gaussian_prior_view(proper, 0, &pv);
    auto selection_status = cosmo_gaussian_selection_view(proper, 0, &history);
    fail_after = -1;
    CHECK(source_status == COSMO_OK && prior_status == COSMO_OK &&
          selection_status == COSMO_OK && live == before_views);
    CHECK(cosmo_gaussian_evaluate(nullptr, &batch, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    CHECK(cosmo_gaussian_evaluate(g, nullptr, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    CHECK(cosmo_gaussian_evaluate(g, &batch, nullptr, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    CHECK(cosmo_gaussian_evaluate(g, &batch, &policy, nullptr) ==
          COSMO_INVALID_INPUT);
    alignas(double) std::array<unsigned char, 32> misaligned{};
    auto alignment_bad = batch;
    alignment_bad.residuals.data =
        reinterpret_cast<const double *>(misaligned.data() + 1);
    CHECK(cosmo_gaussian_evaluate(g, &alignment_bad, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    alignment_bad = batch;
    alignment_bad.ordered_ids.data =
        reinterpret_cast<const cosmo_bytes *>(misaligned.data() + 1);
    CHECK(cosmo_gaussian_evaluate(g, &alignment_bad, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    auto small = policy;
    small.maximum_batch_elements = 1;
    CHECK(cosmo_gaussian_evaluate(g, &batch, &small, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    small = policy;
    small.maximum_string_bytes = 0;
    CHECK(cosmo_gaussian_evaluate(g, &batch, &small, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    small = policy;
    small.maximum_native_bytes = 0;
    CHECK(cosmo_gaussian_evaluate(g, &batch, &small, &out) ==
              COSMO_INVALID_INPUT && !out);
    // Every malformed descriptor fails before borrowed buffers are touched.
    auto malformed = batch;
    malformed.abi_version++;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_ABI_MISMATCH &&
          !out);
    malformed = batch;
    malformed.struct_size--;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    malformed = batch;
    malformed.reserved = 1;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    malformed = batch;
    malformed.mode = 99;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    malformed = batch;
    malformed.row_count = UINT64_MAX;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    malformed = batch;
    malformed.residuals.byte_length--;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    malformed = batch;
    malformed.residuals.element_type = 1;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    malformed = batch;
    malformed.residuals.data = nullptr;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    malformed = batch;
    malformed.ordered_ids.byte_length--;
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) ==
              COSMO_INVALID_INPUT &&
          !out);
    malformed = batch;
    malformed.row_count = 0;
    malformed.residuals = doubles(nullptr, 0);
    CHECK(cosmo_gaussian_evaluate(g, &malformed, &policy, &out) == COSMO_OK &&
          out);
    CHECK(cosmo_gaussian_result_view(out, &rows, &n) == COSMO_OK && n == 0);
    CHECK(cosmo_gaussian_result_destroy(out) == COSMO_OK);
    out = nullptr;
    // Exhaust allocation failure points: owning results/handles never leak and
    // no exception unwinds C.
    bool success = false;
    int points = 0;
    for (int point = 0; point < 200; ++point) {
      auto before = live;
      fail_after = point;
      auto status =
          cosmo_gaussian_proper_offset(g, &prior, &policy, &badg, &semantic, &numerical);
      fail_after = -1;
      if (status == COSMO_OK) {
        CHECK(badg);
        CHECK(cosmo_gaussian_destroy(badg) == COSMO_OK);
        badg = nullptr;
        CHECK(live == before);
        success = true;
        break;
      }
      CHECK(status == COSMO_ALLOCATION_FAILURE && !badg && live == before);
      ++points;
    }
    CHECK(success && points > 0);
    success = false;
    for (int point = 0; point < 100; ++point) {
      auto before = live;
      fail_after = point;
      auto status = cosmo_gaussian_evaluate(g, &batch, &policy, &out);
      fail_after = -1;
      if (status == COSMO_OK) {
        CHECK(out);
        CHECK(cosmo_gaussian_result_destroy(out) == COSMO_OK);
        out = nullptr;
        CHECK(live == before);
        success = true;
        break;
      }
      CHECK(status == COSMO_ALLOCATION_FAILURE && !out && live == before);
    }
    CHECK(success);
    CHECK(cosmo_gaussian_destroy(proper) == COSMO_OK);
    CHECK(cosmo_gaussian_destroy(g) == COSMO_OK);
    CHECK(cosmo_gaussian_destroy(nullptr) == COSMO_OK);
    CHECK(cosmo_gaussian_result_destroy(nullptr) == COSMO_OK);
    CHECK(cosmo_gaussian_source_view(nullptr, &view) == COSMO_INVALID_INPUT);
    CHECK(cosmo_gaussian_result_view(nullptr, &rows, &n) ==
              COSMO_INVALID_INPUT &&
          rows == nullptr && n == 0);
    std::printf("{\"suite\":\"Gaussian_interface_only\",\"checks\":%d,"
                "\"passed\":true}\n",
                checks);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
