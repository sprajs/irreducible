// Independent boundary/science peer for retained piecewise SN. References use
// split composite Simpson in z and an exact rational 3x3 adjugate; production
// analytic segments/ProfileOperator and owner GL8 are not called by the oracle.
#include "irred/piecewise_background.hpp"
#include "irred/supernova.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
namespace {
int checks = 0;
long double maximum_error = 0;
void check(bool b, const char *reason) {
  ++checks;
  if (!b)
    throw std::runtime_error(reason);
}
void close(long double a, long double b) {
  maximum_error = std::max(maximum_error, std::abs(a - b));
  check(std::abs(a - b) <= 1e-10L,
        "frozen synthetic assembled comparison budget");
}
long double expansion(long double z, const std::array<double, 5> &q) {
  long double e = 1;
  const auto &edges = irred::cosmology::piecewise_q_edges;
  for (std::size_t i = 0; i < 5; ++i) {
    const long double lo = edges[i],
                      hi = std::min(z, (long double)edges[i + 1]);
    if (hi <= lo)
      break;
    e *= std::pow((1 + hi) / (1 + lo), 1 + (long double)q[i]);
  }
  return e;
}
long double radial(double z, const std::array<double, 5> &q, int panels) {
  long double total = 0;
  const auto &edges = irred::cosmology::piecewise_q_edges;
  for (std::size_t i = 0; i < 5; ++i) {
    const long double lo = edges[i],
                      hi = std::min((long double)z, (long double)edges[i + 1]);
    if (hi <= lo)
      break;
    const auto h = (hi - lo) / panels;
    long double sum = 0;
    for (int j = 0; j <= panels; ++j)
      sum += (j == 0 || j == panels ? 1
              : j % 2               ? 4
                                    : 2) /
             expansion(lo + j * h, q);
    total += h * sum / 3;
  }
  return total;
}
// C={4,1,0;1,9,2;0,2,16}, det544; exact cofactor inverse.
constexpr long double adj[3][3] = {{140, -16, 2}, {-16, 64, -8}, {2, -8, 35}};
struct ProfileOracle {
  long double coefficient, quadratic;
  std::array<long double, 3> adjusted;
};
ProfileOracle profile(const std::array<long double, 3> &r) {
  const auto coefficient = (126 * r[0] + 40 * r[1] + 29 * r[2]) / 195;
  ProfileOracle out{coefficient, 0, {}};
  for (int i = 0; i < 3; ++i)
    out.adjusted[i] = r[i] - coefficient;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out.quadratic += out.adjusted[i] * adj[i][j] * out.adjusted[j] / 544;
  return out;
}
} // namespace
// Main peer cases are added once the owner's concrete retained-entry API
// freezes.
int main() {
  using namespace irred;
  try {
    observations::Input input;
    input.profile = observations::Profile::gaussian_fixture_v1;
    input.role = observations::Role::synthetic_control;
    input.unit = observations::Unit::magnitude;
    input.calibration = observations::Calibration::not_applicable;
    input.uncertainty = observations::Uncertainty::covariance;
    input.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
    input.component = observations::Component::total;
    input.table_sha256 = std::string(64, 'e');
    input.uncertainty_sha256 = std::string(64, 'f');
    input.ordering_provenance = "independent synthetic correlated coordinates";
    input.measurement_ids = {"one", "two", "three", "excluded"};
    input.event_ids = input.measurement_ids;
    input.uncertainty_axis_ids = input.measurement_ids;
    input.values = {12, 15, 18, 99};
    input.uncertainty_matrix = {4, 1, 0,  8, 1, 9, 2, 7,
                                0, 2, 16, 6, 1, 2, 3, -1};
    input.missing = {0, 0, 0, 0};
    input.quality = {0, 0, 0, 0};
    input.source_selection = {1, 1, 1, 1};
    auto observed = observations::prepare(input, {4, 16, 4096});
    std::array<cosmology::Query, 4> queries{
        {{.1, .11, cosmology::Convention::released_zhd_zhel},
         {.7, .71, cosmology::Convention::released_zhd_zhel},
         {2.5, 2.2, cosmology::Convention::released_zhd_zhel},
         {.01, .01, cosmology::Convention::released_zhd_zhel}}};
    supernova::Policy prep;
    prep.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
    prep.maximum_forward_sensitivity = 1e-10;
    prep.maximum_matrix_elements = 16;
    auto consumer = supernova::prepare_synthetic(
        observed, input.measurement_ids, queries, prep);
    check(consumer.status() == supernova::Status::ok,
          "prepared synthetic correlated selected covariance");
    check(consumer.selected_source_indices().size() == 3,
          "source selection unchanged");
    supernova::PiecewiseEvaluationPolicy policy;
    policy.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
    policy.maximum_forward_sensitivity = 1e-10;
    policy.maximum_models = 6;
    policy.background = {4, 100};
    std::array<supernova::PiecewiseModelPoint, 6> models{
        {supernova::PiecewiseModelPoint({0, 0, 0, 0, 0}),
         supernova::PiecewiseModelPoint({-1, -1, -1, -1, -1}),
         supernova::PiecewiseModelPoint({.5, .5, .5, .5, .5}),
         supernova::PiecewiseModelPoint({-.4, -.4, -.2, .1, .3}),
         supernova::PiecewiseModelPoint({-1, 0, -1, 0, -1}),
         supernova::PiecewiseModelPoint({-3, 2, -3, 2, -3})}};
    auto result = consumer.evaluate_piecewise_batch(models, policy);
    check(result.status == supernova::Status::ok && result.slots.size() == 6,
          "batch owns six ordered rows");
    for (std::size_t m = 0; m < 6; ++m) {
      const auto &s = result.slots[m];
      check(s.status == supernova::Status::ok && s.source.q == models[m].q,
            "source and finite row");
      std::array<long double, 3> residual;
      for (std::size_t i = 0; i < 3; ++i) {
        auto fine = radial(queries[i].z_expansion, models[m].q, 4096),
             coarse = radial(queries[i].z_expansion, models[m].q, 2048);
        check(std::abs(fine - coarse) < 2e-13, "split reference refinement");
        auto prediction =
            5 * std::log10((1 + (long double)queries[i].z_observer) * fine);
        close(s.shape_magnitudes[i], prediction);
        residual[i] = input.values[i] - prediction;
      }
      auto exact = profile(residual);
      close(s.offset_coefficient, exact.coefficient);
      close(s.quadratic, exact.quadratic);
      close(s.relative_profile_score, -exact.quadratic / 2);
      for (int i = 0; i < 3; ++i) {
        close(s.profiled_residuals[i], exact.adjusted[i]);
        check(s.profiled_residuals[i] ==
                  s.solve_diagnostics.adjusted_residuals[i],
              "canonical scored residual retained");
      }
      if (m < 3) {
        supernova::ModelPoint old{cosmology::Model::constant_q_flat_v1, 0,
                                  models[m].q[0]};
        auto same = consumer.evaluate_batch(std::span(&old, 1), prep).slots[0];
        close(s.relative_profile_score, same.relative_profile_score);
      }
    }
    std::size_t total_visits = 0;
    for (const auto &row : result.slots)
      total_visits += row.segment_visits;
    check(result.segment_visits == total_visits,
          "batch work equals all successful row work");
    auto unset_precision = policy;
    unset_precision.arithmetic =
        supernova::PiecewiseEvaluationPolicy{}.arithmetic;
    check(consumer.evaluate_piecewise_batch(models, unset_precision)
              .slots.empty(),
          "new policy does not silently choose legacy precision");
    unset_precision.arithmetic = static_cast<numerics::Arithmetic>(99);
    check(consumer.evaluate_piecewise_batch(models, unset_precision)
              .slots.empty(),
          "unknown precision rejects");
    auto query_cap = policy;
    query_cap.background.maximum_queries = 0;
    auto query_failed =
        consumer
            .evaluate_piecewise_batch(std::span(models.data(), 1), query_cap)
            .slots[0];
    std::printf("query_count_cap background_status %u numerical_status %u\n",
                static_cast<unsigned>(query_failed.background_status),
                static_cast<unsigned>(query_failed.numerical_status));
    check(query_failed.background_status == cosmology::Status::work_limit &&
              query_failed.numerical_status == numerics::Status::work_limit &&
              query_failed.shape_magnitudes.empty(),
          "query-count cap preserves exact work limit cause");
    auto legacy_cap = prep;
    legacy_cap.background.maximum_queries = 0;
    supernova::ModelPoint legacy_point{cosmology::Model::constant_q_flat_v1, 0,
                                       0};
    auto legacy_failure =
        consumer.evaluate_batch(std::span(&legacy_point, 1), legacy_cap)
            .slots[0];
    supernova::ModelPointV2 v2_point{cosmology::Model::constant_q_flat_v1, 0, 0,
                                     -1, 0};
    auto v2_failure =
        consumer.evaluate_batch_v2(std::span(&v2_point, 1), legacy_cap)
            .slots[0]
            .calculation;
    // Read-only witness of inherited cause classification, not a new assertion
    // approving or locking that diagnostic. Legacy production remains
    // unchanged.
    std::printf("legacy_query_cap v1 background %u numerical %u v2 background "
                "%u numerical %u\n",
                static_cast<unsigned>(legacy_failure.background_status),
                static_cast<unsigned>(legacy_failure.numerical_status),
                static_cast<unsigned>(v2_failure.background_status),
                static_cast<unsigned>(v2_failure.numerical_status));
    auto cap = policy;
    cap.background.maximum_segment_visits = result.slots[0].segment_visits;
    auto bounded =
        consumer.evaluate_piecewise_batch(std::span(models.data(), 2), cap);
    check(bounded.slots[0].status == supernova::Status::ok &&
              bounded.slots[1].status != supernova::Status::ok,
          "global segments not reset per model");
    check(bounded.slots[1].shape_magnitudes.empty() &&
              bounded.slots[1].profiled_residuals.empty(),
          "work failure no reusable scientific payload");
    cap.background.maximum_segment_visits = result.slots[0].segment_visits + 1;
    auto partial_work =
        consumer.evaluate_piecewise_batch(std::span(models.data(), 2), cap);
    check(partial_work.slots[1].status != supernova::Status::ok &&
              partial_work.slots[1].segment_visits == 1 &&
              partial_work.segment_visits ==
                  cap.background.maximum_segment_visits &&
              partial_work.segment_visits ==
                  partial_work.slots[0].segment_visits +
                      partial_work.slots[1].segment_visits,
          "batch work includes failed model partial query work");
    for (auto values :
         {std::array<double, 5>{0, 0, -3.0001, 0, 0},
          std::array<double, 5>{0, 0, std::numeric_limits<double>::quiet_NaN(),
                                0, 0}}) {
      supernova::PiecewiseModelPoint bad(values);
      auto row = consumer.evaluate_piecewise_batch(std::span(&bad, 1), policy)
                     .slots[0];
      check(row.status != supernova::Status::ok &&
                row.shape_magnitudes.empty() && row.base_residuals.empty() &&
                row.profiled_residuals.empty(),
            "invalid point fails without arrays");
      check(row.source.q[0] == values[0] && row.source.q[4] == values[4],
            "attempted source endpoints retained");
    }
    auto tight = policy;
    tight.maximum_forward_sensitivity = 1e-30;
    auto failed =
        consumer.evaluate_piecewise_batch(std::span(models.data(), 1), tight)
            .slots[0];
    check(failed.numerical_status ==
                  numerics::Status::conditioning_budget_exceeded &&
              failed.profiled_residuals.empty(),
          "tighter cached policy cause retained");
    cap = policy;
    cap.maximum_models = 0;
    check(consumer.evaluate_piecewise_batch(models, cap).slots.empty(),
          "model cap before row allocation");
    auto beyond = queries;
    beyond[2].z_expansion = std::nextafter(2.5, 3.);
    auto unsupported = supernova::prepare_synthetic(
        observed, input.measurement_ids, beyond, prep);
    auto rejected =
        unsupported
            .evaluate_piecewise_batch(std::span(models.data(), 1), policy)
            .slots[0];
    check(rejected.status != supernova::Status::ok &&
              rejected.shape_magnitudes.empty(),
          "above provider domain cannot drop source row");
    supernova::Consumer unset;
    check(unset.evaluate_piecewise_batch(models, policy).slots.empty(),
          "default consumer no result slots");
    std::printf("piecewise SN independent checks %d PASS max_error %.18Lg\n",
                checks, maximum_error);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL %s after %d\n", e.what(), checks);
    return 1;
  }
}
