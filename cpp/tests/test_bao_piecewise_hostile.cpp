// Independent peer: redshift-segment composite Simpson and exact 2x2
// cofactor Gaussian, versus production analytic segments and Cholesky.
// Observable budget 2e-12+2e-10|reference|; refinement <=10% budget.
// Normalized synthetic density budget1e-8, original24 separately assessed.
#include "irred/bao.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
unsigned checks = 0;
long double worst = 0, refinement = 0, density_error = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
void close(double x, long double y) {
  auto budget = 2e-12L + 2e-10L * std::abs(y);
  worst = std::max(worst, std::abs(x - y) / budget);
  check(std::abs(x - y) <= budget, "observable budget");
}
constexpr std::array<long double, 6> edges{0, .1L, .3L, .6L, 1, 2.5L};
long double expansion(long double z, const std::array<double, 5> &q) {
  long double e = 1;
  for (unsigned k = 0; k < 5; ++k) {
    auto end = std::min(z, edges[k + 1]);
    e *= std::pow((1 + end) / (1 + edges[k]), 1 + (long double)q[k]);
    if (z <= edges[k + 1])
      return e;
  }
  throw std::runtime_error("reference domain");
}
long double integral(double z, const std::array<double, 5> &q, unsigned n) {
  long double sum = 0;
  for (unsigned k = 0; k < 5 && z > edges[k]; ++k) {
    auto lo = edges[k], hi = std::min((long double)z, edges[k + 1]),
         h = (hi - lo) / n, s = 0.L;
    for (unsigned j = 0; j <= n; ++j)
      s += (j == 0 || j == n ? 1 : j % 2 ? 4 : 2) / expansion(lo + j * h, q);
    sum += s * h / 3;
  }
  return sum;
}
long double exact(double z, double q) {
  auto u = 1 + (long double)z;
  if (q == -3)
    return (u * u * u - 1) / 3;
  if (q == 2)
    return (1 - 1 / (u * u)) / 2;
  if (q == 0)
    return std::log(u);
  return (1 - std::pow(u, -q)) / q;
}
std::array<std::array<double, 5>, 8> grid{{{0, 0, 0, 0, 0},
                                           {-1, -1, -1, -1, -1},
                                           {.5, .5, .5, .5, .5},
                                           {-.4, -.4, -.2, .1, .3},
                                           {-1, 0, -1, 0, -1},
                                           {-3, 2, -3, 2, -3},
                                           {-3, -3, -3, -3, -3},
                                           {2, 2, 2, 2, 2}}};
bao::DensityInput input() {
  bao::DensityInput x;
  x.queries = {{.51, bao::Observable::transverse_over_ruler},
               {2.33, bao::Observable::hubble_over_ruler}};
  x.observed = {20, 9};
  x.covariance = {4, 1, 1, 9};
  x.ordered_ids = {"DM-control", "DH-control"};
  x.role = bao::RowRole::synthetic_control;
  x.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  x.table_identity = "independent generated two-coordinate control";
  x.covariance_identity = "integer covariance determinant35";
  x.ordering_provenance = "explicit generated axis order";
  x.calibration_provenance = "unknown";
  x.dependence_provenance = "unknown";
  return x;
}
bao::PreparationPolicy preparation() {
  bao::PreparationPolicy p;
  p.maximum_queries = 2;
  p.maximum_matrix_elements = 4;
  p.maximum_string_bytes = 10000;
  p.maximum_native_bytes = 1000000;
  p.maximum_forward_sensitivity = 1e-10;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  return p;
}
bao::DensityPolicy evaluation() {
  bao::DensityPolicy p;
  p.maximum_models = 30;
  p.maximum_native_bytes = 1000000;
  p.maximum_forward_sensitivity = 1e-10;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.requested = 7;
  p.observables.maximum_queries = 2;
  p.observables.maximum_native_bytes = 1000000;
  p.observables.background.maximum_queries = 2;
  p.observables.background.maximum_native_bytes = 1000000;
  p.observables.background.maximum_segment_visits = 2000;
  return p;
}
} // namespace
int main() {
  try {
    for (auto q : grid)
      for (double ruler : {5000., 10000., 15000.}) {
        auto bg = cosmology::prepare(cosmology::FixedFiveBinQ(q),
                                     cosmology::FlatFLRW{});
        check(bg.status() == cosmology::Status::ok, "prepare reference grid");
        for (double z : {0., .0001, .1, .3, .6, 1., 2.33, 2.5}) {
          std::array<bao::Query, 3> queries{
              {{z, bao::Observable::transverse_over_ruler},
               {z, bao::Observable::hubble_over_ruler},
               {z, bao::Observable::volume_over_ruler}}};
          bao::Policy p;
          p.maximum_queries = 3;
          p.maximum_native_bytes = 1000000;
          p.background.maximum_queries = 3;
          p.background.maximum_native_bytes = 1000000;
          p.background.maximum_segment_visits = 100;
          auto b = bao::evaluate(bg, bao::Ruler(ruler), queries, p);
          check(b.status == cosmology::Status::ok && b.slots.size() == 3,
                "observable batch");
          auto a = integral(z, q, 4096), r = integral(z, q, 8192),
               budget = 2e-12L + 2e-10L * std::abs(r);
          refinement = std::max(refinement, std::abs(a - r) / budget);
          check(std::abs(a - r) <= .1L * budget, "redshift Simpson refinement");
          if (std::all_of(q.begin(), q.end(),
                          [&](double x) { return x == q[0]; })) {
            auto analytic = exact(z, q[0]);
            check(std::abs(r - analytic) <= .1L * budget,
                  "analytic integral control");
            r = analytic;
          }
          const long double scale = 299792.458L / ruler, dm = scale * r,
                            dh = scale / expansion(z, q),
                            dv = std::cbrt(z * dm * dm * dh);
          for (auto &s : b.slots)
            check(s.status == cosmology::Status::ok, "observable finite");
          close(b.slots[0].value.value(), dm);
          close(b.slots[1].value.value(), dh);
          close(b.slots[2].value.value(), dv);
        }
      }
    auto source = input();
    auto edge_bg = cosmology::prepare(cosmology::FixedFiveBinQ(grid[5]),
                                      cosmology::FlatFLRW{});
    std::array<bao::Query, 3> malformed_queries{
        {{.3, static_cast<bao::Observable>(99)},
         {std::nextafter(2.5, std::numeric_limits<double>::infinity()),
          bao::Observable::transverse_over_ruler},
         {.3, bao::Observable::transverse_over_ruler}}};
    bao::Policy edge_policy;
    edge_policy.maximum_queries = 3;
    edge_policy.maximum_native_bytes = 1000000;
    edge_policy.background.maximum_queries = 3;
    edge_policy.background.maximum_native_bytes = 1000000;
    edge_policy.background.maximum_segment_visits = 100;
    auto edge_result = bao::evaluate(edge_bg, bao::Ruler(10000),
                                     malformed_queries, edge_policy);
    check(edge_result.status == cosmology::Status::ok &&
              edge_result.slots[0].status == cosmology::Status::invalid_input &&
              edge_result.slots[1].status ==
                  cosmology::Status::unsupported_domain &&
              edge_result.slots[2].status == cosmology::Status::ok,
          "unknown tag and outside-domain source do not hide valid jump "
          "geometry");
    edge_policy.maximum_queries = 0;
    check(bao::evaluate(edge_bg, bao::Ruler(10000), malformed_queries,
                        edge_policy)
              .slots.empty(),
          "query cap rejects before copies");
    auto dp = preparation();
    auto over_reserved = source;
    over_reserved.covariance.reserve(250000);
    auto oversized = bao::prepare_density(std::move(over_reserved), dp);
    check(oversized.status() != statistics::DensityStatus::finite &&
              oversized.numerical_status() == numerics::Status::work_limit,
          "transferred covariance capacity must fit retained native byte cap");
    for (unsigned kind : {0u, 1u, 2u}) {
      auto capacity_input = source;
      if (kind == 0)
        capacity_input.queries.reserve(100000);
      if (kind == 1)
        capacity_input.ordered_ids[0].reserve(2 * 1024 * 1024);
      if (kind == 2)
        capacity_input.ordering_provenance.reserve(2 * 1024 * 1024);
      auto payload = bao::retained_source_payload_bound(capacity_input);
      check(
          payload && *payload > dp.maximum_native_bytes,
          "reserved native source payload is charged independently of length");
      auto capacity_owner = bao::prepare_density(std::move(capacity_input), dp);
      check(capacity_owner.status() != statistics::DensityStatus::finite &&
                capacity_owner.numerical_status() ==
                    numerics::Status::work_limit,
            "reserved vector or string source cannot bypass retained byte cap");
    }
    auto prepared = bao::prepare_density(source, dp);
    check(prepared.status() == statistics::DensityStatus::finite,
          "density prepared");
    auto policy = evaluation();
    std::vector<bao::ModelPoint> points;
    for (auto q : grid)
      for (double ruler : {5000., 10000., 15000.})
        points.emplace_back(cosmology::FixedFiveBinQ(q), cosmology::FlatFLRW{},
                            bao::Ruler(ruler));
    auto indefinite_source = source;
    indefinite_source.covariance[0] = -4;
    auto diagnostic_owner = bao::prepare_density(indefinite_source, dp);
    check(diagnostic_owner.status() != statistics::DensityStatus::finite &&
              diagnostic_owner.numerical_status() != numerics::Status::ok,
          "indefinite admitted covariance is a numerical preparation failure");
    check(diagnostic_owner.source().ordered_ids == source.ordered_ids &&
              diagnostic_owner.source().table_identity ==
                  source.table_identity &&
              diagnostic_owner.source().covariance ==
                  indefinite_source.covariance,
          "failed factor already retains actual full source identity and "
          "covariance");
    auto retained_failed = diagnostic_owner.retained_payload_bound();
    check(retained_failed && *retained_failed > sizeof(diagnostic_owner),
          "failed diagnostic owner retained source is charged");
    auto diagnostic_result = diagnostic_owner.evaluate(points, policy);
    check(diagnostic_result.status != statistics::DensityStatus::finite &&
              diagnostic_result.slots.empty() &&
              diagnostic_result.numerical_status ==
                  diagnostic_owner.numerical_status(),
          "failed preparation cannot expose a usable density");
    auto out = prepared.evaluate(points, policy);
    check(out.slots.size() == 24, "all24 materialized");
    size_t total = 0;
    for (size_t i = 0; i < 24; ++i) {
      auto &s = out.slots[i];
      check(s.result->density.status == statistics::DensityStatus::finite,
            "density finite");
      check(std::get<cosmology::FixedFiveBinQ>(s.source.expansion).q ==
                    std::get<cosmology::FixedFiveBinQ>(points[i].expansion).q &&
                s.source.ruler.h0_rd_km_s == points[i].ruler.h0_rd_km_s,
            "attempt retained");
      auto q = std::get<cosmology::FixedFiveBinQ>(points[i].expansion).q;
      long double scale = 299792.458L / points[i].ruler.h0_rd_km_s,
                  r0 = 20 - scale * integral(.51, q, 16384),
                  r1 = 9 - scale / expansion(2.33, q);
      auto quad = (9 * r0 * r0 - 2 * r0 * r1 + 4 * r1 * r1) / 35;
      auto target =
          -.5L * (quad + std::log(35.L) + 2 * std::log(2 * std::acos(-1.L)));
      density_error = std::max(density_error,
                               std::abs(s.result->density.log_value - target));
      check(std::abs(s.result->density.log_value - target) <= 1e-8L,
            "cofactor normalized density");
      check(std::abs(s.result->log_determinant - std::log(35.L)) < 1e-13L,
            "determinant normalization");
      check(s.predictions.size() == 2 && s.residuals.size() == 2,
            "owned payload");
      total += s.work.segment_visits;
    }
    check(out.work.segment_visits == total && total <= 2000,
          "global work count");
    auto permuted = source;
    std::swap(permuted.queries[0], permuted.queries[1]);
    std::swap(permuted.observed[0], permuted.observed[1]);
    std::swap(permuted.ordered_ids[0], permuted.ordered_ids[1]);
    permuted.covariance = {9, 1, 1, 4};
    auto permutation_owner = bao::prepare_density(permuted, dp);
    check(permutation_owner.source().ordered_ids == permuted.ordered_ids,
          "matched coordinate permutation preserves declared output order");
    auto permutation_result = permutation_owner.evaluate(points, policy);
    check(permutation_result.slots.size() == out.slots.size(),
          "permutation retains every admitted model");
    for (size_t i = 0; i < out.slots.size(); ++i) {
      check(std::abs(permutation_result.slots[i].result->density.log_value -
                     out.slots[i].result->density.log_value) < 1e-12,
            "full correlated coordinate permutation normalized invariant");
      for (size_t j = 0; j < 2; ++j) {
        check(permutation_result.slots[i].predictions[1-j] == out.slots[i].predictions[j],
              "inverse permutation restores prediction bits");
        check(permutation_result.slots[i].residuals[1-j] == out.slots[i].residuals[j],
              "inverse permutation restores residual bits");
      }
    }
    auto limited = policy;
    limited.observables.background.maximum_segment_visits = 3;
    auto partial = prepared.evaluate(points, limited);
    check(partial.work.segment_visits == 3 &&
              partial.slots[0].result->density.status ==
                  statistics::DensityStatus::finite,
          "first model consumes shared analytic budget");
    check(partial.slots[1].work.segment_visits == 0 &&
              partial.slots[1].numerical_status ==
                  numerics::Status::work_limit &&
              partial.slots[1].predictions.empty(),
          "later model atomic rejection does not reset budget");
    auto tiny = source;
    tiny.queries[0].z = std::numeric_limits<double>::denorm_min();
    auto tiny_owner = bao::prepare_density(tiny, dp);
    check(tiny_owner.status() == statistics::DensityStatus::finite,
          "tiny source structurally retained");
    limited.observables.background.maximum_segment_visits = 1;
    auto failed_work = tiny_owner.evaluate(points, limited);
    check(failed_work.work.segment_visits == 1 &&
              failed_work.slots[0].work.segment_visits == 1,
          "failed computed geometry consumes analytic work");
    check((!failed_work.slots[0].result ||
           failed_work.slots[0].result->density.status !=
               statistics::DensityStatus::finite) &&
              failed_work.slots[0].predictions.empty() &&
              failed_work.slots[0].residuals.empty(),
          "failed geometry cannot publish density payload");
    auto prediction_source = source;
    prediction_source.queries[0].z = 1e-306;
    auto expansion_only =
        cosmology::prepare(points[0].expansion, points[0].geometry);
    auto raw_prediction =
        bao::evaluate(expansion_only, points[0].ruler,
                      prediction_source.queries, policy.observables);
    check(raw_prediction.slots[0].value.has_value(),
          "normal tiny prediction admitted");
    prediction_source.observed[0] =
        std::nextafter(*raw_prediction.slots[0].value,
                       std::numeric_limits<double>::infinity());
    auto prediction_owner = bao::prepare_density(prediction_source, dp);
    auto predictions_only = policy;
    predictions_only.requested = (uint32_t)bao::Output::predictions;
    auto prediction_batch =
        prediction_owner.evaluate(std::span(points).first(1), predictions_only);
    check(prediction_batch.slots[0].predictions.size() == 2 &&
              prediction_batch.slots[0].residuals.empty() &&
              !prediction_batch.slots[0].result,
          "unrequested subnormal residual cannot erase valid predictions");
    check(prediction_batch.slots[0].predictions_state.availability ==
                  cosmology::Availability::available &&
              prediction_batch.slots[0].residuals_state.availability ==
                  cosmology::Availability::not_requested &&
              prediction_batch.slots[0].density_state.availability ==
                  cosmology::Availability::not_requested,
          "predictions-only requested-state truth");
    for (uint32_t mask : {6u, 7u}) {
      predictions_only.requested = mask;
      auto partial_outputs = prediction_owner.evaluate(
          std::span(points).first(1), predictions_only);
      const auto &partial_slot = partial_outputs.slots[0];
      check(partial_slot.predictions.size() == 2 &&
                partial_slot.residuals.empty(),
            "available predictions survive requested residual failure");
      check(partial_slot.predictions_state.availability ==
                    cosmology::Availability::available &&
                partial_slot.residuals_state.availability ==
                    cosmology::Availability::failed &&
                partial_slot.residuals_state.numerical_status ==
                    numerics::Status::outside_domain,
            "residual underflow cause isolated from predictions");
      if (mask == 7) {
        check(
            partial_slot.density_state.availability ==
                    cosmology::Availability::failed &&
                partial_slot.density_state.numerical_status ==
                    numerics::Status::outside_domain &&
                (!partial_slot.result || partial_slot.result->density.status !=
                                             statistics::DensityStatus::finite),
            "density dependency retains underflow failure rather than finite "
            "payload");
      } else {
        check(partial_slot.density_state.availability ==
                      cosmology::Availability::not_requested &&
                  !partial_slot.result,
              "unrequested density remains absent");
      }
    }
    auto bad = policy;
    bad.maximum_models = 0;
    auto capped = prepared.evaluate(points, bad);
    check(capped.slots.empty() &&
              capped.numerical_status == numerics::Status::work_limit,
          "model cap preflight");
    bad = policy;
    bad.maximum_native_bytes = 0;
    capped = prepared.evaluate(points, bad);
    check(capped.slots.empty() &&
              capped.numerical_status == numerics::Status::work_limit,
          "byte cap preflight");
    bao::DensityPolicy unset;
    check(prepared.evaluate(points, unset).slots.empty(),
          "required precision policy");
    bad = policy;
    bad.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
    check(prepared.evaluate(points, bad).slots.empty(),
          "retained precision mismatch");
    auto attempts = points;
    std::get<cosmology::FixedFiveBinQ>(attempts[0].expansion).q[0] = 3;
    attempts[1].ruler.h0_rd_km_s = 4999;
    auto failures = prepared.evaluate(attempts, policy);
    for (size_t i = 0; i < 2; ++i) {
      check(!failures.slots[i].result ||
                failures.slots[i].result->density.status !=
                    statistics::DensityStatus::finite,
            "domain failure");
      check(failures.slots[i].predictions.empty() &&
                failures.slots[i].residuals.empty(),
            "failed payload omitted");
    }
    check(std::get<cosmology::FixedFiveBinQ>(failures.slots[0].source.expansion)
                  .q[0] == 3,
          "invalid source retained");
    auto moved = std::move(prepared);
    check(moved.evaluate(points, policy).slots[0].result->density.status ==
              statistics::DensityStatus::finite,
          "moved owner");
    auto abandoned = prepared.evaluate(points, policy);
    for (const auto &slot : abandoned.slots) {
      check(!slot.result || slot.result->density.status !=
                                statistics::DensityStatus::finite,
            "moved-from cannot publish finite density");
      check(slot.numerical_status != numerics::Status::ok,
            "moved-from failure has explicit numerical cause");
      check(slot.predictions.empty() && slot.residuals.empty(),
            "moved-from cannot publish finite arrays");
    }
    bao::PreparedDensity empty;
    check(empty.evaluate(points, policy).slots.empty(), "default owner");
    std::printf("piecewise BAO peer PASS %u observable_fraction %.18Lg "
                "refinement_fraction %.18Lg density_max_error %.18Lg\n",
                checks, worst, refinement, density_error);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "piecewise BAO peer: %s\n", e.what());
    return 1;
  }
}
