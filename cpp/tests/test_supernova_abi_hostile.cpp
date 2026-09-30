#include "../src/observation_internal.hpp"
#include "irred/supernova.hpp"
#include <bit>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
namespace s = irred::supernova;
namespace o = irred::observations;
namespace c = irred::cosmology;
namespace n = irred::numerics;
namespace {
unsigned checks = 0;
long countdown = -1, live = 0;
void check(bool b, const char *m) {
  ++checks;
  if (!b)
    throw std::runtime_error(m);
}
irred_bytes bytes(const std::string &x) {
  return {reinterpret_cast<const uint8_t *>(x.data()), x.size()};
}
irred_f64_buffer f64(const double *x, size_t z) {
  return {sizeof(irred_f64_buffer), IRRED_ABI_VERSION, 2, 0, x, z,
          z * sizeof(double)};
}
struct Fixture {
  std::string
      a = "a",
      b = std::string(100, 'b'), hash = std::string(64, 'c'),
      chash = std::string(64, 'd'),
      order = "explicit synthetic numeric control on measured magnitude scale";
  irred_bytes ids[2]{bytes(a), bytes(b)};
  double values[2]{20, 21}, C[4]{4, 1, 1, 9}, q = -1, epsilon = .2;
  uint8_t missing[2]{0, 0};
  uint64_t quality[2]{0, 0}, indices[2]{0, 1};
  irred_observation_descriptor source{};
  irred_observation_policy observation_policy{2, 4, 4096};
  irred_magnitude_coordinate coords[2]{
      {sizeof(irred_magnitude_coordinate), IRRED_ABI_VERSION, .1, .12, 1, 0},
      {sizeof(irred_magnitude_coordinate), IRRED_ABI_VERSION, .2, .25, 1, 0}};
  irred_magnitude_selection selection{};
  irred_supernova_preparation_policy prep{};
  irred_supernova_evaluation_policy eval{};
  irred_supernova_model model{};
  irred_supernova_batch batch{};
  Fixture() {
    source.struct_size = sizeof(source);
    source.abi_version = IRRED_ABI_VERSION;
    source.profile = (uint32_t)o::Profile::typed_magnitude_covariance;
    source.role = (uint32_t)o::Role::observed_measurement;
    source.unit = (uint32_t)o::Unit::magnitude;
    source.calibration = (uint32_t)o::Calibration::unknown;
    source.uncertainty = (uint32_t)o::Uncertainty::covariance;
    source.uncertainty_unit = (uint32_t)o::UncertaintyUnit::magnitude_squared;
    source.component = (uint32_t)o::Component::total;
    source.table_sha256 = bytes(hash);
    source.uncertainty_sha256 = bytes(chash);
    source.ordering_provenance = bytes(order);
    source.measurement_ids = source.event_ids =
        source.uncertainty_axis_ids = {ids, 2, sizeof(ids)};
    source.values = f64(values, 2);
    source.uncertainty_matrix = f64(C, 4);
    source.zhd = source.zcmb = source.zhel = f64(nullptr, 0);
    source.missing = {missing, 2, 2};
    source.quality = {quality, 2, sizeof(quality)};
    selection = {sizeof(selection),
                 IRRED_ABI_VERSION,
                 0,
                 0,
                 {indices, 2, sizeof(indices)},
                 coords,
                 2,
                 sizeof(coords)};
    prep = {sizeof(prep), IRRED_ABI_VERSION, 1, 0, 2, 4, 4096, 1000000, 1e-10};
    eval.struct_size = sizeof(eval);
    eval.abi_version = IRRED_ABI_VERSION;
    eval.arithmetic = 1;
    eval.requested = 63;
    eval.maximum_models = 2;
    eval.maximum_array_elements = 100;
    eval.maximum_native_bytes = 1000000;
    eval.maximum_forward_sensitivity = 1e-10;
    eval.projection = {sizeof(irred_projection_policy),
                       IRRED_ABI_VERSION,
                       1,
                       24,
                       1e-12,
                       1e-12,
                       100000,
                       2,
                       100000,
                       100};
    model = {
        sizeof(model),
        IRRED_ABI_VERSION,
        0,
        0,
        {sizeof(irred_expansion_spec), IRRED_ABI_VERSION, 1, 0, f64(&q, 1)},
        {sizeof(irred_source_effect_spec), IRRED_ABI_VERSION, 1, 0,
         f64(&epsilon, 1)}};
    batch = {sizeof(batch), IRRED_ABI_VERSION, &model, 1, sizeof(model)};
  }
};
struct View {
  const irred_supernova_row *rows = nullptr;
  uint64_t count = 0, callbacks = 0, segments = 0;
  uint32_t status = 0, num = 0;
  void get(irred_supernova_result *r) {
    check(irred_supernova_result_view(r, &rows, &count, &status, &num,
                                              &callbacks,
                                              &segments) == IRRED_OK,
          "view");
  }
};
} // namespace
void *operator new(size_t z) {
  if (countdown == 0)
    throw std::bad_alloc();
  if (countdown > 0)
    --countdown;
  auto p = std::malloc(z ? z : 1);
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
void operator delete(void *p, size_t) noexcept { operator delete(p); }
int main() {
  try {
    {
      Fixture bad;
      bad.C[0] = -4;
      irred_prepared *raw = nullptr;
      uint32_t raw_status = 0;
      check(irred_prepare_observations(&bad.source, &bad.observation_policy,
                                       &raw, &raw_status) == IRRED_OK &&
                raw_status == 0,
            "indefinite finite source structurally admitted");
      irred_supernova *failed = nullptr;
      check(irred_supernova_prepare(raw, &bad.selection, &bad.prep,
                                            &failed) == IRRED_OK &&
                failed,
            "failed factor owns admitted diagnostic source");
      irred_supernova_view failed_view{};
      check(irred_supernova_source_view(failed, &failed_view) ==
                    IRRED_OK &&
                failed_view.status != 0,
            "failed factor source view reports numerical failure");
      check(failed_view.ordered_ids.length == 2 &&
                failed_view.selected_source_indices.length == 2 &&
                failed_view.coordinate_count == 2,
            "failed factor retains exact attempted selected subset");
      check(
          failed_view.selected_source_indices.data[0] == 0 &&
              failed_view.selected_source_indices.data[1] == 1 &&
              failed_view.coordinates[0].z_expansion == .1 &&
              failed_view.coordinates[1].observer_redshift == .25,
          "failed selected source coordinates and order remain interpretable");
      bad.eval.maximum_native_bytes = 0;
      irred_supernova_result *diagnostic = nullptr;
      check(irred_supernova_evaluate(failed, &bad.batch, &bad.eval,
                                             &diagnostic) == IRRED_OK &&
                diagnostic,
            "failed prepared cause takes precedence over configured evaluation "
            "quota");
      View failed_result;
      failed_result.get(diagnostic);
      check(failed_result.count == 0 && failed_result.status != 0 &&
                failed_result.num == failed_view.preparation_numerical_status,
            "failed owner cannot expose a usable profile or false quota cause");
      irred_supernova_destroy(failed);
      irred_observation_destroy(raw);
      uint32_t has_source = 0;
      check(irred_supernova_result_source_view(
                diagnostic, &failed_view, &has_source) == IRRED_OK &&
                has_source && failed_view.coordinate_count == 2 &&
                failed_view.source.uncertainty_matrix.data[0] == -4,
            "failed result retains raw and selected source after all parent "
            "handles release");
      irred_supernova_result_destroy(diagnostic);
    }
    Fixture f;
    irred_prepared *obs = nullptr;
    uint32_t semantic = 0;
    check(irred_prepare_observations(&f.source, &f.observation_policy, &obs,
                                     &semantic) == IRRED_OK &&
              semantic == 0,
          "measured observation");
    const long prep_baseline = live;
    unsigned prep_failures = 0;
    for (long i = 0; i < 200; ++i) {
      irred_supernova *probe = nullptr;
      countdown = i;
      const auto rc =
          irred_supernova_prepare(obs, &f.selection, &f.prep, &probe);
      countdown = -1;
      if (rc == IRRED_OK) {
        check(probe != nullptr, "prepare allocation sweep success owner");
        irred_supernova_destroy(probe);
        check(live == prep_baseline, "prepare success cleanup");
        break;
      }
      check(rc == IRRED_ALLOCATION_FAILURE && !probe,
            "prepare allocation failure null");
      check(live == prep_baseline, "prepare failure cleanup");
      ++prep_failures;
    }
    check(prep_failures > 5, "prepare allocation coverage");
    irred_supernova *owner = nullptr;
    check(irred_supernova_prepare(obs, &f.selection, &f.prep, &owner) ==
                  IRRED_OK &&
              owner,
          "generic prepare");
    irred_supernova_view source{};
    check(irred_supernova_source_view(owner, &source) == IRRED_OK &&
              source.status == 0 && source.source.role == 0,
          "measured source preserved");
    auto shared = shared_native_observations(obs);
    s::SelectedMagnitudeSource selected{
        shared,
        {0, 1},
        {f.a, f.b},
        {{.1, {.12, c::Convention::released_zhd_zhel}},
         {.2, {.25, c::Convention::released_zhd_zhel}}}};
    auto native = s::prepare(
        selected, {n::Arithmetic::longdouble_cpu_v1, 2, 4, 1e-10, 1000000});
    s::Policy p;
    p.arithmetic = n::Arithmetic::longdouble_cpu_v1;
    p.requested = 63;
    p.maximum_models = 2;
    p.maximum_native_bytes = 1000000;
    p.maximum_forward_sensitivity = 1e-10;
    p.background.maximum_queries = 2;
    p.background.maximum_native_bytes = 1000000;
    p.background.maximum_callbacks = 100000;
    p.background.maximum_segment_visits = 100;
    p.background.integration = n::IntegrationPolicy{1e-12, 1e-12, 100000, 24};
    s::ModelPoint point(c::ConstantQ(-1), c::FlatFLRW{},
                        s::GreyLog1pMagnitude(.2));
    auto direct = native.evaluate_batch(std::span(&point, 1), p);
    irred_supernova_result *result = nullptr;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
                  IRRED_OK &&
              result,
          "evaluate");
    View v;
    v.get(result);
    check(v.count == 1 && v.rows[0].score.state.availability == 1,
          "profile available");
    check(std::bit_cast<uint64_t>(v.rows[0].score.relative_profile_score) ==
              std::bit_cast<uint64_t>(
                  direct.slots[0].score->relative_profile_score),
          "native score bits");
    for (size_t i = 0; i < 2; ++i) {
      check(std::bit_cast<uint64_t>(v.rows[0].geometric_shape.data[i]) ==
                std::bit_cast<uint64_t>(direct.slots[0].geometric_shape[i]),
            "geometry bits");
      check(std::bit_cast<uint64_t>(v.rows[0].magnitude_effect.data[i]) ==
                std::bit_cast<uint64_t>(direct.slots[0].magnitude_shifts[i]),
            "effect bits");
    }
    irred_observation_destroy(obs);
    obs = nullptr;
    irred_supernova_destroy(owner);
    owner = nullptr;
    f.epsilon = -.2;
    f.values[0] = 99;
    f.b.assign(100, 'z');
    uint32_t available = 0;
    check(irred_supernova_result_source_view(result, &source,
                                                     &available) == IRRED_OK &&
              available,
          "result source after parent destruction");
    check(source.source.values.data[0] == 20 &&
              source.source.measurement_ids.data[1].data[0] == 'b',
          "owned source lifetime");
    check(v.rows[0].source.source_effect.parameters.data[0] == .2,
          "attempted effect lifetime");
    irred_supernova_result_destroy(result);
    check(irred_prepare_observations(&f.source, &f.observation_policy, &obs,
                                     &semantic) == IRRED_OK &&
              semantic == 0,
          "second observations");
    check(irred_supernova_prepare(obs, &f.selection, &f.prep, &owner) ==
                  IRRED_OK &&
              owner,
          "second prepare");
    f.eval.maximum_array_elements = 8;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
              IRRED_OK,
          "four F64 vector groups exact array bound");
    v.get(result);
    check(v.count == 1 && v.rows[0].score.state.availability == 1,
          "fixed diagnostics and U64 node mappings are not extra F64 vectors");
    irred_supernova_result_destroy(result);
    f.eval.maximum_array_elements = 7;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
              IRRED_OK,
          "array bound minus one diagnostic owner");
    v.get(result);
    check(v.count == 0 && v.num == (uint32_t)n::Status::work_limit,
          "array quota rejects entire result before compute");
    irred_supernova_result_destroy(result);
    f.eval.maximum_array_elements = 100;
    f.eval.projection.maximum_callbacks = 0;
    f.epsilon = std::numeric_limits<double>::denorm_min();
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
              IRRED_OK,
          "simultaneous independent scientific failures owned");
    v.get(result);
    check(v.count == 1 && v.rows[0].geometry_state.availability == 3 &&
              v.rows[0].geometry_state.numerical_status ==
                  (uint32_t)n::Status::work_limit,
          "geometry retains its own work limit cause");
    check(v.rows[0].effect_state.availability == 3 &&
              v.rows[0].effect_state.numerical_status ==
                  (uint32_t)n::Status::outside_domain,
          "effect retains its own underflow cause");
    check(v.rows[0].corrected_residuals.length == 0 &&
              v.rows[0].score.state.availability == 3,
          "dependent failed results have no usable payload");
    irred_supernova_result_destroy(result);
    f.epsilon = -.2;
    f.eval.requested = 4;
    f.eval.projection.has_integration = 0;
    f.eval.projection.maximum_depth = 0;
    f.eval.projection.absolute_tolerance =
        f.eval.projection.relative_tolerance = 0;
    f.eval.projection.maximum_evaluations_per_integral =
        f.eval.projection.maximum_callbacks =
            f.eval.projection.maximum_segment_visits = 0;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
              IRRED_OK,
          "B only");
    v.get(result);
    check(v.count == 1 && v.rows[0].effect_state.availability == 1 &&
              v.callbacks == 0 && v.segments == 0,
          "B only no quadrature");
    check(v.rows[0].geometry_state.availability == 0, "geometry not requested");
    irred_supernova_result_destroy(result);
    f.q = 3;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
              IRRED_OK,
          "invalid model scientific owner");
    v.get(result);
    check(v.count == 1 && v.rows[0].effect_state.availability == 3 &&
              v.rows[0].magnitude_effect.length == 0,
          "B only still admits full model domain");
    irred_supernova_result_destroy(result);
    f.q = -1;
    f.eval.maximum_native_bytes = 0;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
                  IRRED_OK &&
              result,
          "quota owner");
    v.get(result);
    check(v.count == 0 && v.num == (uint32_t)n::Status::work_limit,
          "empty quota");
    irred_supernova_result_destroy(result);
    f.eval.maximum_native_bytes = 1000000;
    f.eval.maximum_native_bytes = 0;
    f.eval.arithmetic = 0;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
              IRRED_OK,
          "valid precision mismatch semantic owner before quota");
    v.get(result);
    check(v.count == 0 &&
              v.status == (uint32_t)s::Status::incompatible_metadata &&
              v.num != (uint32_t)n::Status::work_limit,
          "precision mismatch not a configured resource failure");
    check(irred_supernova_result_source_view(result, &source,
                                                     &available) == IRRED_OK &&
              available && source.arithmetic == 1,
          "mismatch exposes actual prepared arithmetic");
    irred_supernova_result_destroy(result);
    f.eval.arithmetic = 2;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
                  IRRED_INVALID_INPUT &&
              !result,
          "unknown precision is structural null");
    f.eval.arithmetic = 1;
    f.eval.maximum_native_bytes = 1000000;
    f.batch.model_count = UINT64_MAX;
    f.batch.model_byte_length = 0;
    check(irred_supernova_evaluate(owner, &f.batch, &f.eval, &result) ==
                  IRRED_INVALID_INPUT &&
              !result,
          "overflow descriptor null");
    f.batch.model_count = 1;
    f.batch.model_byte_length = sizeof(f.model);
    long before = live;
    unsigned failures = 0;
    for (long i = 0; i < 200; ++i) {
      countdown = i;
      auto rc =
          irred_supernova_evaluate(owner, &f.batch, &f.eval, &result);
      countdown = -1;
      if (rc == IRRED_OK) {
        check(result != nullptr, "allocation success");
        irred_supernova_result_destroy(result);
        check(live == before, "success cleanup");
        break;
      }
      check(rc == IRRED_ALLOCATION_FAILURE && !result, "allocation null");
      check(live == before, "failure cleanup");
      ++failures;
    }
    check(failures > 5, "allocation coverage");
    irred_supernova_destroy(owner);
    irred_observation_destroy(obs);
    {
      Fixture hard;
      std::vector<std::string> names(17);
      std::vector<irred_bytes> ids;
      std::vector<double> values(17, 20), covariance(17 * 17, 0);
      std::vector<uint8_t> missing(17, 0);
      std::vector<uint64_t> quality(17, 0), indices;
      std::vector<irred_magnitude_coordinate> coordinates;
      for (size_t i = 0; i < 17; ++i) {
        names[i] = "hard-row-" + std::to_string(i);
        ids.push_back(bytes(names[i]));
        covariance[i * 17 + i] = 1;
        indices.push_back(i);
        coordinates.push_back(hard.coords[0]);
      }
      hard.source.measurement_ids = hard.source.event_ids =
          hard.source.uncertainty_axis_ids = {ids.data(), 17, 17 * sizeof(irred_bytes)};
      hard.source.values = f64(values.data(), 17);
      hard.source.uncertainty_matrix = f64(covariance.data(), 289);
      hard.source.missing = {missing.data(), 17, 17};
      hard.source.quality = {quality.data(), 17, 17 * sizeof(uint64_t)};
      hard.observation_policy = {17, 289, 4096};
      uint32_t hard_status = 0;
      check(irred_prepare_observations(&hard.source, &hard.observation_policy,
                                      &obs, &hard_status) == IRRED_OK && hard_status == 0, "hard-count source");
      std::vector<irred_supernova_model> models(65536, hard.model);
      hard.batch = {sizeof(hard.batch), IRRED_ABI_VERSION, models.data(),
                    models.size(), models.size() * sizeof(models[0])};
      hard.prep.maximum_selected_rows = 17;
      hard.prep.maximum_matrix_elements = 289;
      hard.eval.maximum_native_bytes = 0;
      hard.eval.requested = 1;
      for (size_t selected : {16u, 17u}) {
        hard.selection.source_indices = {indices.data(), selected, selected * sizeof(uint64_t)};
        hard.selection.coordinates = coordinates.data();
        hard.selection.coordinate_count = selected;
        hard.selection.coordinate_byte_length = selected * sizeof(coordinates[0]);
        check(irred_supernova_prepare(obs, &hard.selection, &hard.prep,
                                              &owner) == IRRED_OK,
              "hard-count consumer preparation");
        const auto rc = irred_supernova_evaluate(owner, &hard.batch,
                                                        &hard.eval, &result);
        if (selected == 16) {
          check(rc == IRRED_OK && result, "exact hard M*N admitted before configured quota");
          v.get(result);
          check(v.count == 0 && v.num == (uint32_t)n::Status::work_limit,
                "hard boundary diagnostic has zero result count");
          irred_supernova_result_destroy(result);
        } else
          check(rc == IRRED_INVALID_INPUT && !result,
                "hard M*N exceeded rejects before result allocation");
        irred_supernova_destroy(owner);
      }
      hard.batch.models = nullptr;
      check(irred_supernova_evaluate(nullptr, &hard.batch, &hard.eval,
                                            &result) == IRRED_INVALID_INPUT && !result,
            "null descriptors reset result without dereference");
      irred_observation_destroy(obs);
    }
    std::cout << "Current SN ABI peer PASS " << checks
              << " allocation failures " << failures
              << " preparation allocation failures " << prep_failures << '\n';
  } catch (const std::exception &x) {
    countdown = -1;
    std::cerr << "FAIL " << checks << ": " << x.what() << '\n';
    return 1;
  }
}
