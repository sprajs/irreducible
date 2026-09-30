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
  const std::array<supernova::MagnitudeCoordinate, 2> coordinates{
      {{1, {.9, cosmology::Convention::released_zhd_zhel}},
       {2, {2.1, cosmology::Convention::released_zhd_zhel}}}};
  auto selected = [&](const observations::Prepared &value) {
    auto owner = std::make_shared<const observations::Prepared>(value);
    return supernova::SelectedMagnitudeSource{
        owner,
        {0, 1},
        {"synthetic:A", "synthetic:B"},
        {coordinates.begin(), coordinates.end()}};
  };
  supernova::PreparationPolicy preparation;
  preparation.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  preparation.maximum_selected_rows = 3;
  preparation.maximum_matrix_elements = 9;
  preparation.maximum_forward_sensitivity = 1e-10;
  preparation.maximum_native_bytes = 512 * 1024 * 1024;
  supernova::Policy policy;
  policy.arithmetic = preparation.arithmetic;
  policy.maximum_forward_sensitivity = 1e-10;
  policy.maximum_models = 8;
  policy.maximum_native_bytes = 1 << 20;
  policy.requested = 63;
  policy.background.integration =
      numerics::IntegrationPolicy{1e-14, 1e-14, 100000, 30};
  policy.background.maximum_queries = 3;
  policy.background.maximum_callbacks = 200000;
  policy.background.maximum_segment_visits = 100;
  policy.background.maximum_native_bytes = 1 << 20;
  const auto peak = supernova::preparation_payload_bound(
      original, 2, preparation.arithmetic);
  check(peak.has_value(), "representable preparation peak");
  auto bounded = preparation;
  bounded.maximum_native_bytes = *peak;
  check(supernova::prepare(selected(original), bounded).status() ==
            supernova::Status::ok,
        "exact preparation peak admits");
  --bounded.maximum_native_bytes;
  check(supernova::prepare(selected(original), bounded).status() ==
            supernova::Status::work_limit,
        "preparation peak minus one rejects");
  auto reserved = selected(original);
  reserved.source_indices.reserve(100000);
  bounded.maximum_native_bytes = *peak;
  check(supernova::prepare(std::move(reserved), bounded).status() ==
            supernova::Status::work_limit,
        "transferred selected reserve capacity counted");
  check(!supernova::preparation_payload_bound(original, SIZE_MAX,
                                              preparation.arithmetic),
        "preparation overflow rejects");
  auto consumer = supernova::prepare(selected(original), preparation);
  check(consumer.status() == supernova::Status::ok,
        "selected synthetic consumer");
  check(consumer.selected_source().source_indices.size() == 2,
        "same strict threshold");
  std::array<supernova::ModelPoint, 1> models{
      {supernova::ModelPoint{cosmology::ConstantQ{-1}, cosmology::FlatFLRW{},
                             supernova::NoMagnitudeEffect{}}}};
  auto result = consumer.evaluate_batch(models, policy);
  check(result.status == supernova::Status::ok && result.slots.size() == 1 &&
            result.slots[0].status == supernova::Status::ok,
        "analytic model batch");
  const auto &slot = result.slots[0];
  near(slot.score->offset_coefficient, 10.L + 7.L / 11, 1e-10,
       "offset adjugate oracle");
  near(slot.score->quadratic, 25.L / 11, 1e-10, "profile quadratic oracle");
  near(slot.score->relative_profile_score, -25.L / 22, 1e-10,
       "relative score not density");
  near(slot.geometric_shape[0], 5 * ln(1.9L) / ln(10), 1e-10,
       "zHEL prefactor convention");
  check(consumer.selected_source().source->source().uncertainty_matrix ==
            original.source().uncertainty_matrix,
        "excluded raw asymmetry retained");
  const auto &adjusted = slot.profiled_residuals;
  const long double direct_q =
      (9.L * adjusted[0] * adjusted[0] - 2.L * adjusted[0] * adjusted[1] +
       4.L * adjusted[1] * adjusted[1]) /
      35;
  near(slot.score->quadratic, direct_q, 1e-10,
       "reported score uses retained canonical adjusted residuals");
  auto large_input = original.source();
  large_input.values[0] += 1e8;
  large_input.values[1] += 1e8;
  auto large_source =
      observations::prepare(std::move(large_input), {3, 9, 4096});
  auto large_consumer = supernova::prepare(selected(large_source), preparation);
  auto large = large_consumer.evaluate_batch(models, policy);
  check(large.slots[0].status == supernova::Status::ok,
        "large shared intercept remains profile score");
  const auto difference =
      static_cast<long double>(large.slots[0].corrected_residuals[0]) -
      large.slots[0].corrected_residuals[1];
  near(large.slots[0].score->quadratic, difference * difference / 11, 1e-10,
       "large intercept stable difference oracle not subtraction of large "
       "quadratics");
  auto wrong = selected(original);
  std::swap(wrong.ordered_ids[0], wrong.ordered_ids[1]);
  check(supernova::prepare(std::move(wrong), preparation).status() ==
            supernova::Status::incompatible_metadata,
        "explicit supplied order mismatch");
  check(supernova::prepare(selected(original), {}).status() ==
            supernova::Status::invalid_input,
        "untouched preparation policy invalid");
  auto missing_policy = preparation;
  missing_policy.maximum_forward_sensitivity = 0;
  check(supernova::prepare(selected(original), missing_policy).status() ==
            supernova::Status::invalid_input,
        "explicit solver budget required");
  auto cap = policy;
  cap.maximum_models = 0;
  check(consumer.evaluate_batch(models, cap).slots.empty(),
        "batch cap before slot allocation");
  std::array<supernova::ModelPoint, 2> repeated{{models[0], models[0]}};
  auto work = policy;
  work.background.maximum_callbacks = 1;
  auto limited = consumer.evaluate_batch(repeated, work);
  check(limited.slots.size() == 2 &&
            limited.slots[0].numerical_status == numerics::Status::work_limit &&
            limited.slots[1].numerical_status == numerics::Status::work_limit,
        "whole batch work cap truthful failure");
  auto global = policy;
  global.background.maximum_callbacks = slot.work.callbacks + 2;
  auto limited_after_first = consumer.evaluate_batch(repeated, global);
  check(limited_after_first.slots[0].status == supernova::Status::ok &&
            limited_after_first.slots[1].numerical_status ==
                numerics::Status::work_limit,
        "callback allowance is not reset per model");
  check(limited_after_first.slots[0].work.callbacks +
                limited_after_first.slots[1].work.callbacks <=
            global.background.maximum_callbacks,
        "batch callback accounting within cap");
  auto expansion =
      cosmology::prepare(cosmology::ConstantQ{-1}, cosmology::FlatFLRW{});
  std::array<cosmology::Request, 1> request{
      {cosmology::Request{1, (uint32_t)cosmology::Observable::luminosity_shape,
                          coordinates[0].observer}}};
  auto shape = expansion.evaluate(request, policy.background);
  near(5 * std::log10(*shape.slots[0].luminosity_shape.value),
       slot.geometric_shape[0], 1e-10, "shape requires no computational H0");
  std::printf("{\"suite\":\"supernova_owner\",\"checks\":%d,\"passed\":true}\n",
              checks);
}
