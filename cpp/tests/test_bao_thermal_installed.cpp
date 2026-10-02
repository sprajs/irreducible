// Standalone consumer candidate: compiled separately against an installed
// archive and headers as well as the source-tree CTest target.
#include <iostream>
#include <irred/bao_thermal.hpp>
int main() {
  using namespace irred;
  bao::DensityInput input;
  input.queries = {{0, bao::Observable::transverse_over_ruler}};
  input.observed = {0};
  input.covariance = {1};
  input.ordered_ids = {"synthetic-DM0"};
  input.role = bao::RowRole::synthetic_control;
  input.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  input.table_identity = "synthetic zero-distance control";
  input.covariance_identity = "unit variance";
  input.ordering_provenance = "single explicit row";
  input.calibration_provenance = "supplied drag";
  input.dependence_provenance = "single synthetic source";
  auto observations = bao::prepare_density(
      std::move(input),
      {1, 1, 4096, 1024 * 1024, 1e-8, numerics::Arithmetic::longdouble_cpu_v1});
  const cosmology::ThermalObservableRequest request{
      {70, .0245, .1225, 2.7255, 1e-5, {}},
      1059.95,
      "synthetic supplied drag",
      "explicit physical density source"};
  bao::ThermalDensityPolicy policy;
  policy.maximum_models = 1;
  policy.maximum_queries = 1;
  policy.maximum_string_bytes = 4096;
  policy.maximum_native_bytes = 1024 * 1024;
  policy.maximum_total_callbacks = 100000000;
  policy.maximum_forward_sensitivity = 1e-8;
  policy.requested = 7;
  const auto result =
      observations.evaluate_thermal(std::span(&request, 1), policy);
  if (result.slots.size() != 1 || !result.slots[0].result ||
      result.slots[0].predictions.size() != 1 ||
      result.slots[0].predictions[0] != 0 ||
      result.slots[0].result->density.status !=
          statistics::DensityStatus::finite ||
      result.callbacks != result.outer_callbacks + result.momentum_callbacks)
    return 1;
  std::cout << "PASS installed thermal BAO consumer\n";
}
