#include "irred/abi.h"
#include "irred/bao.hpp"
#include <bit>
#include <cstdlib>
#include <iostream>
#include <new>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace b = irred::bao;
namespace c = irred::cosmology;
namespace n = irred::numerics;
namespace {
unsigned checks = 0;
long countdown = -1, live = 0;
void check(bool x, const char *why) {
  ++checks;
  if (!x)
    throw std::runtime_error(why);
}
cosmo_bytes bytes(const std::string &s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
cosmo_f64_buffer f64(const double *p, size_t z) {
  return {sizeof(cosmo_f64_buffer), COSMO_ABI_VERSION, 2, 0, p, z,
          z * sizeof(double)};
}
struct Fixture {
  std::string a = "short", longid = std::string(100, 'b'),
              table = std::string(64, 'a'), covhash = std::string(64, 'c'),
              order = "synthetic explicit axis order",
              cal = "no calibration asserted",
              dep = "declared correlated control, no joint target",
              redshift = "P01/released-effective-redshift/v1",
              ruler = std::string(b::ruler_convention_id);
  cosmo_bytes ids[2]{bytes(a), bytes(longid)};
  cosmo_bao_query queries[2]{{.1, 0, 0}, {.2, 1, 0}};
  double y[2]{4, 32}, C[4]{4, 1, 1, 9}, q = -1;
  cosmo_bao_source src{};
  cosmo_current_bao_preparation_policy prep{};
  cosmo_current_bao_evaluation_policy eval{};
  cosmo_current_bao_model model{};
  cosmo_current_bao_batch batch{};
  Fixture() {
    src = {sizeof(src),
           COSMO_ABI_VERSION,
           1,
           0,
           queries,
           2,
           sizeof(queries),
           f64(y, 2),
           f64(C, 4),
           {ids, 2, sizeof(ids)},
           {ids, 2, sizeof(ids)},
           bytes(table),
           bytes(covhash),
           bytes(order),
           bytes(cal),
           bytes(dep),
           bytes(redshift),
           bytes(ruler)};
    prep = {sizeof(prep), COSMO_ABI_VERSION, 1, 0, 2, 4, 4096, 1000000, 1e-10};
    eval.struct_size = sizeof(eval);
    eval.abi_version = COSMO_ABI_VERSION;
    eval.arithmetic = 1;
    eval.requested = 7;
    eval.projection.struct_size = sizeof(eval.projection);
    eval.projection.abi_version = COSMO_ABI_VERSION;
    eval.projection.has_integration = 1;
    eval.projection.maximum_depth = 24;
    eval.projection.absolute_tolerance = 1e-12;
    eval.projection.relative_tolerance = 1e-12;
    eval.projection.maximum_evaluations_per_integral = 100000;
    eval.projection.maximum_queries = 2;
    eval.projection.maximum_callbacks = 100000;
    eval.projection.maximum_segment_visits = 100;
    eval.maximum_models = 2;
    eval.maximum_array_elements = 100;
    eval.maximum_native_bytes = 1000000;
    eval.maximum_forward_sensitivity = 1e-10;
    model = {
        sizeof(model),
        COSMO_ABI_VERSION,
        0,
        0,
        {sizeof(cosmo_expansion_spec), COSMO_ABI_VERSION, 1, 0, f64(&q, 1)},
        10000};
    batch = {sizeof(batch), COSMO_ABI_VERSION, &model, 1, sizeof(model)};
  }
  b::DensityInput native() {
    b::DensityInput x;
    x.queries = {{.1, b::Observable::transverse_over_ruler},
                 {.2, b::Observable::hubble_over_ruler}};
    x.observed = {y[0], y[1]};
    x.covariance = {4, 1, 1, 9};
    x.ordered_ids = {a, longid};
    x.role = b::RowRole::synthetic_control;
    x.covariance_unit = b::CovarianceUnit::dimensionless_ratio_squared;
    x.table_identity = table;
    x.covariance_identity = covhash;
    x.ordering_provenance = order;
    x.calibration_provenance = cal;
    x.dependence_provenance = dep;
    return x;
  }
};
struct View {
  const cosmo_current_bao_row *rows = nullptr;
  uint64_t count = 0, callbacks = 0, segments = 0;
  uint32_t status = 0, num = 0;
  void get(cosmo_current_bao_result *r) {
    check(cosmo_current_bao_result_view(r, &rows, &count, &status, &num,
                                        &callbacks, &segments) == COSMO_OK,
          "result view");
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
    Fixture f;
    cosmo_current_bao *owner = nullptr;
    cosmo_current_bao_result *result = nullptr;
    f.queries[0].z = -1;
    const auto rejected = cosmo_current_bao_prepare(&f.src, &f.prep, &owner);
    check(rejected == COSMO_INVALID_INPUT && !owner,
          "invalid source domain is structural rejection without empty owner");
    f.queries[0].z = .1;
    f.C[0] = -4;
    check(cosmo_current_bao_prepare(&f.src, &f.prep, &owner) == COSMO_OK && owner,
          "admitted indefinite covariance retains diagnostic owner");
    cosmo_current_bao_view failed_source{};
    check(cosmo_current_bao_source_view(owner, &failed_source) == COSMO_OK &&
              failed_source.source.query_count == 2 &&
              failed_source.source.covariance.data[0] == -4,
          "failed factor retains exact nonempty source");
    cosmo_current_bao_destroy(owner);
    owner = nullptr;
    f.C[0] = 4;
    cosmo_bytes reversed_ids[2]{f.ids[1], f.ids[0]};
    f.src.ordered_ids = {reversed_ids, 2, sizeof(reversed_ids)};
    check(cosmo_current_bao_prepare(&f.src, &f.prep, &owner) ==
                  COSMO_INVALID_INPUT && !owner,
          "row-only permutation rejects covariance axis mismatch");
    f.src.ordered_ids = {f.ids, 2, sizeof(f.ids)};
    f.src.covariance_axis_ids = {reversed_ids, 2, sizeof(reversed_ids)};
    check(cosmo_current_bao_prepare(&f.src, &f.prep, &owner) ==
                  COSMO_INVALID_INPUT && !owner,
          "axis-only permutation rejects covariance axis mismatch");
    f.src.covariance_axis_ids = {f.ids, 2, sizeof(f.ids)};
    check(cosmo_current_bao_prepare(&f.src, &f.prep, &owner) == COSMO_OK &&
              owner,
          "prepare");
    cosmo_current_bao_view source{};
    check(cosmo_current_bao_source_view(owner, &source) == COSMO_OK &&
              source.status == 0,
          "source view");
    auto native =
        b::prepare_density(f.native(), {2, 4, 4096, 1000000, 1e-10,
                                        n::Arithmetic::longdouble_cpu_v1});
    b::DensityPolicy p;
    p.arithmetic = n::Arithmetic::longdouble_cpu_v1;
    p.maximum_forward_sensitivity = 1e-10;
    p.maximum_models = 2;
    p.maximum_native_bytes = 1000000;
    p.requested = 7;
    p.observables.maximum_queries = 2;
    p.observables.maximum_native_bytes = 1000000;
    p.observables.background.maximum_queries = 2;
    p.observables.background.maximum_native_bytes = 1000000;
    p.observables.background.maximum_callbacks = 100000;
    p.observables.background.maximum_segment_visits = 100;
    p.observables.background.integration =
        n::IntegrationPolicy{1e-12, 1e-12, 100000, 24};
    b::ModelPoint point(c::ConstantQ(-1), c::FlatFLRW{}, b::Ruler(10000));
    auto direct = native.evaluate(std::span(&point, 1), p);
    check(cosmo_current_bao_evaluate(owner, &f.batch, &f.eval, &result) ==
                  COSMO_OK &&
              result,
          "evaluate");
    View v;
    v.get(result);
    check(v.count == 1 && v.rows[0].density.state.availability == 1,
          "finite density state");
    check(
        std::bit_cast<uint64_t>(v.rows[0].density.log_density) ==
            std::bit_cast<uint64_t>(direct.slots[0].result->density.log_value),
        "exact native density parity");
    for (size_t i = 0; i < 2; ++i) {
      check(std::bit_cast<uint64_t>(v.rows[0].predictions.data[i]) ==
                std::bit_cast<uint64_t>(direct.slots[0].predictions[i]),
            "prediction bits");
      check(std::bit_cast<uint64_t>(v.rows[0].residuals.data[i]) ==
                std::bit_cast<uint64_t>(direct.slots[0].residuals[i]),
            "residual bits");
    }
    cosmo_current_bao_destroy(owner);
    owner = nullptr;
    f.y[0] = 9;
    f.q = 0;
    f.longid.assign(100, 'z');
    uint32_t available = 0;
    check(cosmo_current_bao_result_source_view(result, &source, &available) ==
                  COSMO_OK &&
              available,
          "source outlives prepare");
    check(source.source.observed.data[0] == 4 &&
              source.source.ordered_ids.data[1].data[0] == 'b',
          "owned source bytes");
    check(v.rows[0].source.expansion.parameters.data[0] == -1,
          "attempted model bits owned");
    cosmo_current_bao_result_destroy(result);
    {
      Fixture partial;
      partial.queries[0].z = 1e-306;
      auto tiny_queries = partial.native().queries;
      tiny_queries[0].z = 1e-306;
      auto expansion = c::prepare(c::ConstantQ(-1), c::FlatFLRW{});
      auto tiny = b::evaluate(expansion, b::Ruler(10000), tiny_queries,
                              p.observables);
      check(tiny.slots[0].value.has_value(), "tiny normal prediction admitted");
      partial.y[0] = std::nextafter(*tiny.slots[0].value,
                                    std::numeric_limits<double>::infinity());
      cosmo_current_bao *partial_owner = nullptr;
      check(cosmo_current_bao_prepare(&partial.src, &partial.prep,
                                      &partial_owner) == COSMO_OK,
            "tiny residual source admitted");
      for (uint32_t mask : {2u, 6u, 7u}) {
        partial.eval.requested = mask;
        check(cosmo_current_bao_evaluate(partial_owner, &partial.batch,
                                         &partial.eval, &result) == COSMO_OK,
              "partial requested outputs materialized");
        View pv;
        pv.get(result);
        check(pv.count == 1 && pv.rows[0].predictions_state.availability == 1 &&
                  pv.rows[0].predictions.length == 2,
              "available prediction survives residual underflow");
        check(pv.rows[0].residuals_state.availability == (mask == 2 ? (uint32_t)c::Availability::not_requested : (uint32_t)c::Availability::failed),
              "unrequested or failed residual truth");
        if (mask != 2)
          check(pv.rows[0].residuals_state.numerical_status ==
                    (uint32_t)n::Status::outside_domain &&
                    pv.rows[0].residuals.length == 0,
                "residual underflow cause and payload omission");
        if (mask == 7)
          check(pv.rows[0].density.state.availability == (uint32_t)c::Availability::failed,
                "failed density not usable despite owned optional result");
        cosmo_current_bao_result_destroy(result);
      }
      cosmo_current_bao_destroy(partial_owner);
    }
    check(cosmo_current_bao_prepare(&f.src, &f.prep, &owner) == COSMO_OK,
          "second preparation");
    f.eval.maximum_native_bytes = 0;
    check(cosmo_current_bao_evaluate(owner, &f.batch, &f.eval, &result) ==
                  COSMO_OK &&
              result,
          "quota diagnostic");
    v.get(result);
    check(v.count == 0 && v.num == (uint32_t)n::Status::work_limit,
          "empty quota worklimit");
    cosmo_current_bao_result_destroy(result);
    f.eval.maximum_native_bytes = 1000000;
    f.batch.model_count = UINT64_MAX;
    f.batch.model_byte_length = 0;
    check(cosmo_current_bao_evaluate(owner, &f.batch, &f.eval, &result) ==
                  COSMO_INVALID_INPUT &&
              !result,
          "impossible descriptor before reading");
    f.batch.model_count = 1;
    f.batch.model_byte_length = sizeof(f.model);
    f.eval.abi_version = 0;
    check(cosmo_current_bao_evaluate(owner, &f.batch, &f.eval, &result) ==
                  COSMO_ABI_MISMATCH &&
              !result,
          "version null reset");
    f.eval.abi_version = COSMO_ABI_VERSION;
    long before = live;
    unsigned failures = 0;
    for (long i = 0; i < 200; ++i) {
      countdown = i;
      auto rc = cosmo_current_bao_evaluate(owner, &f.batch, &f.eval, &result);
      countdown = -1;
      if (rc == COSMO_OK) {
        check(result != nullptr, "sweep success");
        cosmo_current_bao_result_destroy(result);
        check(live == before, "success cleanup");
        break;
      }
      check(rc == COSMO_ALLOCATION_FAILURE && !result, "sweep allocation null");
      check(live == before, "failure cleanup");
      ++failures;
    }
    check(failures > 5, "allocation coverage");
    cosmo_current_bao_destroy(owner);
    owner = nullptr;
    before = live;
    unsigned preparation_failures = 0;
    for (long i = 0; i < 200; ++i) {
      countdown = i;
      const auto rc = cosmo_current_bao_prepare(&f.src, &f.prep, &owner);
      countdown = -1;
      if (rc == COSMO_OK) {
        check(owner != nullptr, "preparation sweep success owner");
        cosmo_current_bao_destroy(owner);
        owner = nullptr;
        check(live == before, "preparation success cleanup");
        break;
      }
      check(rc == COSMO_ALLOCATION_FAILURE && !owner,
            "preparation allocation failure null owner");
      check(live == before, "preparation allocation failure cleanup");
      ++preparation_failures;
    }
    check(preparation_failures > 10, "preparation allocation coverage");
    {
      Fixture hard;
      std::vector<cosmo_current_bao_model> models(65536, hard.model);
      hard.batch = {sizeof(hard.batch), COSMO_ABI_VERSION, models.data(),
                    models.size(), models.size() * sizeof(models[0])};
      hard.eval.requested = 1;
      hard.eval.maximum_native_bytes = 0;
      hard.prep.maximum_queries = 17;
      hard.prep.maximum_matrix_elements = 289;
      for (size_t rows : {16u, 17u}) {
        std::vector<std::string> names(rows);
        std::vector<cosmo_bytes> ids;
        std::vector<cosmo_bao_query> queries(rows, {.1, 0, 0});
        std::vector<double> y(rows, 4), covariance(rows * rows, 0);
        for (size_t i = 0; i < rows; ++i) {
          names[i] = "hard-row-" + std::to_string(i);
          ids.push_back(bytes(names[i]));
          covariance[i * rows + i] = 1;
        }
        hard.src.queries = queries.data();
        hard.src.query_count = rows;
        hard.src.query_byte_length = rows * sizeof(queries[0]);
        hard.src.observed = f64(y.data(), rows);
        hard.src.covariance = f64(covariance.data(), covariance.size());
        hard.src.ordered_ids = hard.src.covariance_axis_ids =
            {ids.data(), rows, rows * sizeof(ids[0])};
        check(cosmo_current_bao_prepare(&hard.src, &hard.prep, &owner) == COSMO_OK,
              "hard-count preparation");
        const auto rc = cosmo_current_bao_evaluate(owner, &hard.batch,
                                                  &hard.eval, &result);
        if (rows == 16) {
          check(rc == COSMO_OK && result,
                "exact hard M*N admitted before configured quota");
          v.get(result);
          check(v.count == 0 && v.num == (uint32_t)n::Status::work_limit,
                "hard boundary empty diagnostic result count");
          cosmo_current_bao_result_destroy(result);
          hard.batch.models = nullptr;
          check(cosmo_current_bao_evaluate(owner, &hard.batch, &hard.eval,
                                           &result) == COSMO_INVALID_INPUT && !result,
                "null model buffer rejects before reading");
          hard.batch.models = models.data();
        } else
          check(rc == COSMO_INVALID_INPUT && !result,
                "hard M*N excess structural null before allocation");
        cosmo_current_bao_destroy(owner);
        owner = nullptr;
      }
    }
    std::cout << "Current BAO ABI peer PASS " << checks
              << " allocation failures " << failures
              << " preparation allocation failures " << preparation_failures
              << '\n';
  } catch (const std::exception &x) {
    countdown = -1;
    std::cerr << "FAIL " << checks << ": " << x.what() << '\n';
    return 1;
  }
}
