// Native synthetic consumer comparison: de Sitter I=z and rational adjugate
// covariance [[4,1],[1,9]]. Profile residual (2,-3) has a=7/11,
// q=25/11. Magnitude logarithms use an independent long-double atanh series.
// Named synthetic absolute score/offset budget1e-10; no real-release claim.
#include "irred/supernova.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
using namespace irred;
namespace {
int checks = 0;
void check(bool ok, const char *why) {
  ++checks;
  if (!ok)
    throw std::runtime_error(why);
}
void near(double actual, long double expected, double budget, const char *why) {
  check(std::abs(actual - expected) <= budget, why);
}
long double ln(long double x) {
  int exponent = 0;
  while (x >= 2) {
    x /= 2;
    ++exponent;
  }
  while (x < 1) {
    x *= 2;
    --exponent;
  }
  auto series = [](long double z) {
    long double sum = 0, p = z;
    for (int k = 1; k < 240; k += 2) {
      sum += p / k;
      p *= z * z;
    }
    return 2 * sum;
  };
  return series((x - 1) / (x + 1)) + exponent * series(1.L / 3);
}
observations::Prepared source() {
  observations::Input s{};
  s.profile = observations::Profile::gaussian_fixture_v1;
  s.role = observations::Role::synthetic_control;
  s.unit = observations::Unit::magnitude;
  s.calibration = observations::Calibration::not_applicable;
  s.uncertainty = observations::Uncertainty::covariance;
  s.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  s.component = observations::Component::total;
  s.ordering_provenance = "explicit generated ordered measurement keys";
  s.table_sha256 = std::string(64, 'a');
  s.uncertainty_sha256 = std::string(64, 'b');
  s.measurement_ids = {"synthetic:A", "synthetic:B", "synthetic:excluded"};
  s.event_ids = {"event:A", "event:B", "event:excluded"};
  s.uncertainty_axis_ids = s.measurement_ids;
  s.values = {static_cast<double>(5 * ln(1.9L) / ln(10) + 12),
              static_cast<double>(5 * ln(6.2L) / ln(10) + 7), 0};
  s.missing = {0, 0, 0};
  s.quality = {0, 0, 0};
  s.source_selection = {1, 1, 1};
  s.uncertainty_matrix = {4, 1, 8, 1, 9, 5, 7, 5, -100};
  return observations::prepare(std::move(s), {3, 9, 4096});
}
} // namespace
int main() {
  auto original = source();
  check(original.status() == observations::Status::ok,
        "truthful synthetic source");
  std::array<cosmology::Query, 3> queries{
      {{1, .9, cosmology::Convention::released_zhd_zhel},
       {2, 2.1, cosmology::Convention::released_zhd_zhel},
       {.001, .001, cosmology::Convention::released_zhd_zhel}}};
  supernova::Policy policy;
  policy.maximum_forward_sensitivity = 1e-10;
  policy.maximum_matrix_elements = 9;
  auto consumer = supernova::prepare_synthetic(
      original, original.source().measurement_ids, queries, policy);
  check(consumer.status() == supernova::Status::ok,
        "selected synthetic consumer");
  check(consumer.selected_source_indices().size() == 2,
        "same strict threshold");
  std::array<supernova::ModelPoint, 1> models{
      {{cosmology::Model::constant_q_flat_v1, 0, -1}}};
  auto result = consumer.evaluate_batch(models, policy);
  check(result.status == supernova::Status::ok && result.slots.size() == 1 &&
            result.slots[0].status == supernova::Status::ok,
        "analytic model batch");
  const auto &slot = result.slots[0];
  near(slot.offset_coefficient, 10.L + 7.L / 11, 1e-10,
       "offset adjugate oracle");
  near(slot.quadratic, 25.L / 11, 1e-10, "profile quadratic oracle");
  near(slot.relative_profile_score, -25.L / 22, 1e-10,
       "relative score not density");
  near(slot.shape_magnitudes[0], 5 * ln(1.9L) / ln(10), 1e-10,
       "zHEL prefactor convention");
  check(consumer.observations().source().uncertainty_matrix ==
            original.source().uncertainty_matrix,
        "excluded raw asymmetry retained");
  const auto &adjusted = slot.profiled_residuals;
  const long double direct_q =
      (9.L * adjusted[0] * adjusted[0] - 2.L * adjusted[0] * adjusted[1] +
       4.L * adjusted[1] * adjusted[1]) /
      35;
  near(slot.quadratic, direct_q, 1e-10,
       "reported score uses retained canonical adjusted residuals");
  auto large_input = original.source();
  large_input.values[0] += 1e8;
  large_input.values[1] += 1e8;
  auto large_source =
      observations::prepare(std::move(large_input), {3, 9, 4096});
  auto large_consumer = supernova::prepare_synthetic(
      large_source, large_source.source().measurement_ids, queries, policy);
  auto large = large_consumer.evaluate_batch(models, policy);
  check(large.slots[0].status == supernova::Status::ok,
        "large shared intercept remains profile score");
  const auto difference =
      static_cast<long double>(large.slots[0].base_residuals[0]) -
      large.slots[0].base_residuals[1];
  near(large.slots[0].quadratic, difference * difference / 11, 1e-10,
       "large intercept stable difference oracle not subtraction of large "
       "quadratics");
  auto wrong = original.source().measurement_ids;
  std::swap(wrong[0], wrong[1]);
  check(
      supernova::prepare_synthetic(original, wrong, queries, policy).status() ==
          supernova::Status::incompatible_metadata,
      "explicit supplied order mismatch");
  check(supernova::prepare_synthetic(
            original, original.source().measurement_ids, queries, {})
                .status() == supernova::Status::invalid_input,
        "untouched policy has no default scientific sensitivity");
  auto missing_policy = policy;
  missing_policy.maximum_forward_sensitivity = 0;
  check(supernova::prepare_synthetic(original,
                                     original.source().measurement_ids, queries,
                                     missing_policy)
                .status() == supernova::Status::invalid_input,
        "explicit solver budget required");
  auto cap = policy;
  cap.maximum_models = 0;
  check(consumer.evaluate_batch(models, cap).slots.empty(),
        "batch cap before slot allocation");
  std::array<supernova::ModelPoint, 2> repeated{{models[0], models[0]}};
  auto work = policy;
  work.background.maximum_total_evaluations = 1;
  auto limited = consumer.evaluate_batch(repeated, work);
  check(limited.slots.size() == 2 &&
            limited.slots[0].status == supernova::Status::work_limit &&
            limited.slots[1].status == supernova::Status::work_limit,
        "whole batch work cap truthful failure");
  auto global = policy;
  global.background.maximum_total_evaluations = slot.background_evaluations + 2;
  auto limited_after_first = consumer.evaluate_batch(repeated, global);
  check(limited_after_first.slots[0].status == supernova::Status::ok &&
            limited_after_first.slots[1].status ==
                supernova::Status::work_limit,
        "callback allowance is not reset per model");
  check(limited_after_first.slots[0].background_evaluations +
                limited_after_first.slots[1].background_evaluations <=
            global.background.maximum_total_evaluations,
        "batch callback accounting within cap");
  for (double h0 : {40., 100.}) {
    auto background =
        cosmology::prepare({cosmology::Model::constant_q_flat_v1, h0, 0, -1});
    auto b = background.evaluate_batch(queries, policy.background);
    near(5 * std::log10(b.slots[0].dimensionless_luminosity_shape),
         slot.shape_magnitudes[0], 1e-10, "intercept free H0 invariance");
  }
  std::printf("{\"suite\":\"supernova_owner\",\"checks\":%d,\"passed\":true}\n",
              checks);
}
