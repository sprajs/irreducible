// Independent consumer contract checks: synthetic inputs only, exact de Sitter
// I=z and diagonal C=diag(2,2). Profiling two residuals gives q=(r0-r1)^2/4.
// Shared libm logs not independent transcendental qualification. Score budget
// 1e-10 predeclared; no original-data or full W01 reproduction claim.
#include "irred/supernova.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
using namespace irred;
int main() {
  try {
    int checks = 0;
    auto check = [&](bool v) {
      ++checks;
      if (!v)
        throw std::runtime_error("consumer hostile contract");
    };
    observations::Input s;
    s.profile = observations::Profile::gaussian_fixture_v1;
    s.role = observations::Role::synthetic_control;
    s.unit = observations::Unit::magnitude;
    s.calibration = observations::Calibration::not_applicable;
    s.uncertainty = observations::Uncertainty::covariance;
    s.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
    s.component = observations::Component::total;
    s.ordering_provenance = "generated ordered axes";
    s.table_sha256 = std::string(64, 'a');
    s.uncertainty_sha256 = std::string(64, 'b');
    s.measurement_ids = {"first", "second", "threshold"};
    s.event_ids = s.measurement_ids;
    s.uncertainty_axis_ids = s.measurement_ids;
    s.values = {3, 7, 999};
    s.missing = {0, 0, 0};
    s.quality = {0, 0, 0};
    s.source_selection = {1, 1, 1};
    s.uncertainty_matrix = {2, 0, 99, 0, 2, 99, -99, -99, -1};
    auto data = observations::prepare(s, {3, 9, 4096});
    check(data.status() == observations::Status::ok);
    std::array<cosmology::Query, 3> q{
        {{1, 1, cosmology::Convention::released_zhd_zhel},
         {2, 2, cosmology::Convention::released_zhd_zhel},
         {.01, .01, cosmology::Convention::released_zhd_zhel}}};
    supernova::Policy p;
    p.maximum_forward_sensitivity = 1e-10;
    auto c = supernova::prepare_synthetic(data, s.measurement_ids, q, p);
    check(c.status() == supernova::Status::ok);
    check(c.selected_source_indices().size() == 2);
    check(c.probability_metadata().matrix_validation_scope ==
          statistics::MatrixValidationScope::selected_covariance_only);
    std::array<supernova::ModelPoint, 3> models{
        {{cosmology::Model::constant_q_flat_v1, 0, -1},
         {static_cast<cosmology::Model>(999), 0, 0},
         {cosmology::Model::constant_q_flat_v1, 0, -1}}};
    auto result = c.evaluate_batch(models, p);
    check(result.slots.size() == 3);
    check(result.slots[0].status == supernova::Status::ok);
    check(result.slots[1].status != supernova::Status::ok);
    check(result.slots[1].shape_magnitudes.empty() &&
          result.slots[1].profiled_residuals.empty());
    check(result.slots[2].status == supernova::Status::ok);
    long double r0 = 3 - 5 * std::log10(2.L), r1 = 7 - 5 * std::log10(6.L);
    check(std::abs(result.slots[0].quadratic - (r0 - r1) * (r0 - r1) / 4) <
          1e-10L);
    auto wrong = q;
    wrong[0].convention = cosmology::Convention::geometric_same_redshift;
    check(supernova::prepare_synthetic(data, s.measurement_ids, wrong, p)
              .status() == supernova::Status::incompatible_metadata);
    auto tighter = p;
    tighter.maximum_forward_sensitivity =
        c.profile_operator().cached_response_forward_sensitivity() / 2;
    auto failed = c.evaluate_batch(
        std::span<const supernova::ModelPoint>(models.data(), 1), tighter);
    check(failed.slots[0].status == supernova::Status::numerical_failure);
    check(failed.slots[0].numerical_status ==
          numerics::Status::conditioning_budget_exceeded);
    check(failed.slots[0].base_residuals.empty() &&
          failed.slots[0].profiled_residuals.empty());
    std::printf("{\"suite\":\"supernova_independent_consumer_contract\","
                "\"checks\":%d,\"passed\":true}\n",
                checks);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
