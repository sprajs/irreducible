// Independent analytic BAO ratio fixtures (flat de Sitter, EdS, constant-q)
// derive from direct antiderivatives; no astronomy/reference runtime
// dependency.
#include "fixtures/bao_reference.hpp"
#include "fixtures/ldlt_reference.hpp"
#include "irred/bao.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace irred;
namespace {
int checks = 0;
double worst = 0;
void check(bool ok, const char *why) {
  ++checks;
  if (!ok)
    throw std::runtime_error(why);
}
void close(double actual, long double expected) {
  const auto difference = std::abs(static_cast<long double>(actual) - expected);
  const auto budget = 2e-12L + 2e-10L * std::abs(expected);
  worst = std::max(worst, static_cast<double>(difference / budget));
  check(difference <= budget, "analytic observable budget");
}
cosmology::Expansion cq(double q) {
  return cosmology::prepare(cosmology::ConstantQ{q}, cosmology::FlatFLRW{});
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 1)
      throw std::runtime_error(
          "original11 moved to test_bao_reference optional guard");
    (void)argv;
    bao::Policy p;
    p.background.integration =
        numerics::IntegrationPolicy{1e-14, 1e-13, 100000, 30};
    p.maximum_queries = 100;
    p.maximum_native_bytes = 1 << 20;
    p.background.maximum_queries = 100;
    p.background.maximum_callbacks = 2000000;
    p.background.maximum_segment_visits = 2000;
    p.background.maximum_native_bytes = 1 << 20;
    for (double q : {-1., -.5, -1e-10, 0., 1e-10, .5})
      for (double z : {0., .0001, .295, .51, 1.484, 2.33, 5.})
        for (double hrd : {5000., 10000., 15000.}) {
          std::vector<bao::Query> queries{
              {z, bao::Observable::transverse_over_ruler},
              {z, bao::Observable::hubble_over_ruler},
              {z, bao::Observable::volume_over_ruler}};
          auto b = bao::evaluate(cq(q), bao::Ruler(hrd), queries, p);
          check(b.status == cosmology::Status::ok && b.slots.size() == 3,
                "analytic batch");
          const auto l = std::log1p(static_cast<long double>(z));
          const auto integral = q == 0 ? l : -std::expm1(-q * l) / q;
          const auto dm = 299792.458L / hrd * integral;
          const auto dh = 299792.458L / hrd * std::exp(-(1 + q) * l);
          const long double expected[]{dm, dh, std::cbrt(z * dm * dm * dh)};
          for (std::size_t i = 0; i < 3; ++i) {
            check(b.slots[i].status == cosmology::Status::ok,
                  "analytic finite");
            check(b.slots[i].numerical_status == numerics::Status::ok,
                  "underlying finite");
            close(*b.slots[i].value, expected[i]);
          }
        }
    std::vector<bao::Query> qs{{.1, bao::Observable::transverse_over_ruler},
                               {1., bao::Observable::hubble_over_ruler},
                               {2.33, bao::Observable::volume_over_ruler}};
    auto ref = bao::evaluate(cq(-.5), bao::Ruler(10000), qs, p);
    // Computational H0 no longer exists in an expansion hypothesis or BAO
    // projection; the free-ruler scaling control below preserves its algebra.
    auto doubled = bao::evaluate(cq(-.5), bao::Ruler(5000), qs, p);
    for (std::size_t i = 0; i < 3; ++i)
      close(*doubled.slots[i].value, 2.L * *ref.slots[i].value);
    auto lambda =
        cosmology::prepare(cosmology::CPL{.3, -1, 0}, cosmology::FlatFLRW{});
    auto lcdm = cosmology::prepare(cosmology::LCDM{.3}, cosmology::FlatFLRW{});
    auto l = bao::evaluate(lambda, bao::Ruler(10000), qs, p);
    auto m = bao::evaluate(lcdm, bao::Ruler(10000), qs, p);
    for (std::size_t i = 0; i < 3; ++i)
      check(l.slots[i].value == m.slots[i].value, "Lambda exact limit");
    qs.push_back({1, static_cast<bao::Observable>(99)});
    qs.push_back({-1, bao::Observable::transverse_over_ruler});
    qs.push_back({std::numeric_limits<double>::quiet_NaN(),
                  bao::Observable::hubble_over_ruler});
    auto bad = bao::evaluate(cq(0), bao::Ruler(10000), qs, p);
    check(bad.slots[3].status == cosmology::Status::invalid_input,
          "unknown tag rejected");
    check(bad.slots[4].status == cosmology::Status::unsupported_domain,
          "negative z rejected");
    check(bad.slots[5].status == cosmology::Status::invalid_input,
          "NaN rejected");
    check(!bad.slots[3].value && !bad.slots[4].value, "no failed payload");
    for (double hrd :
         {0., 4999., 15001., std::numeric_limits<double>::infinity()}) {
      auto b = bao::evaluate(cq(0), bao::Ruler(hrd), qs, p);
      check(b.status != cosmology::Status::ok && b.slots.empty(),
            "ruler invalid empty");
    }
    auto small = p;
    small.maximum_queries = 1;
    check(bao::evaluate(cq(0), bao::Ruler(10000), qs, small).slots.empty(),
          "count cap before copy");
    small = p;
    small.maximum_native_bytes = sizeof(bao::Slot) - 1;
    check(bao::evaluate(cq(0), bao::Ruler(10000), qs, small).status ==
              cosmology::Status::work_limit,
          "byte cap");
    small = p;
    small.background.maximum_callbacks = 3;
    auto capped = bao::evaluate(cq(0), bao::Ruler(10000), qs, small);
    check(capped.work.callbacks <= 3, "global callback cap");
    check(capped.slots[0].status == cosmology::Status::work_limit &&
              capped.slots[1].value &&
              capped.slots[2].status == cosmology::Status::work_limit,
          "DH standalone E survives depleted integration budget");
    auto empty = bao::evaluate(cq(0), bao::Ruler(10000), {}, p);
    check(empty.status == cosmology::Status::ok && empty.slots.empty() &&
              empty.work.callbacks == 0,
          "empty batch");
    bao::DensityInput input;
    input.queries = {{.1, bao::Observable::transverse_over_ruler},
                     {.2, bao::Observable::hubble_over_ruler}};
    input.observed = {299792.458 / 10000 * .1 + 1, 299792.458 / 10000 + 2};
    input.covariance = {4, 1, 1, 9};
    input.ordered_ids = {"synthetic:DM:.1", "synthetic:DH:.2"};
    input.role = bao::RowRole::synthetic_control;
    input.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
    input.table_identity = "analytic deSitter synthetic ratios";
    input.covariance_identity = "rational covariance determinant35";
    input.ordering_provenance = "explicit fixture row order";
    input.calibration_provenance = "synthetic no calibration";
    input.dependence_provenance =
        "two correlated synthetic coordinates; no crossprobe claim";
    bao::PreparationPolicy preparation{
        100, 4, 10000, 100000, 1e-10, numerics::Arithmetic::longdouble_cpu_v1};
    bao::DensityPolicy dp{
        p, 10, 100000, 1e-10, numerics::Arithmetic::longdouble_cpu_v1, 7};
    // Deliberately no braces: public default initialization must be safe.
    bao::DensityInput unset_input;
    bao::PreparationPolicy unset_policy;
    bao::PreparedDensity unset_prepared;
    check(unset_input.role == bao::RowRole::unknown &&
              unset_input.covariance_unit == bao::CovarianceUnit::unknown,
          "unset source metadata is unknown");
    check(unset_policy.maximum_queries == 0 &&
              unset_policy.maximum_matrix_elements == 0 &&
              unset_policy.maximum_string_bytes == 0 &&
              unset_policy.maximum_native_bytes == 0 &&
              unset_policy.maximum_forward_sensitivity == 0 &&
              unset_policy.arithmetic !=
                  numerics::Arithmetic::binary64_legacy_v1 &&
              unset_policy.arithmetic !=
                  numerics::Arithmetic::longdouble_cpu_v1,
          "unset policy cannot silently choose arithmetic or resources");
    check(unset_prepared.source().role == bao::RowRole::unknown &&
              unset_prepared.source().covariance_unit ==
                  bao::CovarianceUnit::unknown,
          "default prepared source safely inspectable");
    auto unset_failure = bao::prepare_density(input, unset_policy);
    check(unset_failure.status() != statistics::DensityStatus::finite &&
              unset_failure.source().role == bao::RowRole::unknown,
          "unset required policy rejects without fabricated released source");
    for (int field = 0; field < 6; ++field) {
      auto incomplete = preparation;
      switch (field) {
      case 0:
        incomplete.maximum_queries = 0;
        break;
      case 1:
        incomplete.maximum_matrix_elements = 0;
        break;
      case 2:
        incomplete.maximum_string_bytes = 0;
        break;
      case 3:
        incomplete.maximum_native_bytes = 0;
        break;
      case 4:
        incomplete.maximum_forward_sensitivity = 0;
        break;
      case 5:
        incomplete.arithmetic = unset_policy.arithmetic;
        break;
      }
      auto rejected = bao::prepare_density(input, incomplete);
      check(rejected.status() != statistics::DensityStatus::finite &&
                rejected.source().role == bao::RowRole::unknown,
            "unusable preparation field rejects before retained source");
    }
    auto partial = input;
    partial.role = unset_input.role;
    check(bao::prepare_density(partial, preparation).status() !=
              statistics::DensityStatus::finite,
          "unset role rejects despite otherwise explicit valid input");
    partial = input;
    partial.covariance_unit = unset_input.covariance_unit;
    check(bao::prepare_density(partial, preparation).status() !=
              statistics::DensityStatus::finite,
          "unset covariance unit rejects");
    auto prepared = bao::prepare_density(input, preparation);
    check(prepared.status() == statistics::DensityStatus::finite,
          "density prepares");
    std::vector<bao::ModelPoint> models{
        {cosmology::ConstantQ{-1}, cosmology::FlatFLRW{}, bao::Ruler(10000)}};
    auto density = prepared.evaluate(models, dp);
    check(density.slots.size() == 1 &&
              density.slots[0].result->density.status ==
                  statistics::DensityStatus::finite,
          "density finite");
    close(density.slots[0].result->quadratic, .6L);
    const auto expected_log =
        -.5L * (.6L + std::log(35.L) + 2 * std::log(2 * std::acos(-1.L)));
    check(std::abs(density.slots[0].result->density.log_value - expected_log) <
              1e-8,
          "normalized determinant density budget");
    check(density.slots[0].predictions.size() == 2 &&
              density.slots[0].residuals.size() == 2,
          "retained ordered arrays");
    auto tightdp = dp;
    tightdp.maximum_forward_sensitivity = 1e-30;
    auto failed = prepared.evaluate(models, tightdp);
    check(failed.slots[0].numerical_status ==
              numerics::Status::conditioning_budget_exceeded,
          "tight exact failure cause");
    check(failed.slots[0].predictions.size() == 2 &&
              failed.slots[0].residuals.size() == 2,
          "requested predictions survive independent density failure");
    models.emplace_back(cosmology::ConstantQ{-1}, cosmology::FlatFLRW{},
                        bao::Ruler(0));
    auto mixed = prepared.evaluate(models, dp);
    check(
        mixed.slots[1].source.ruler.h0_rd_km_s == 0 &&
            std::get<cosmology::ConstantQ>(mixed.slots[1].source.expansion).q ==
                -1,
        "attempted active identity");
    check(mixed.slots[1].result &&
              mixed.slots[1].result->density.status !=
                  statistics::DensityStatus::finite &&
              mixed.slots[1].predictions.empty(),
          "invalid ruler no finite density");
    auto original = input;
    input.ordered_ids[1] = input.ordered_ids[0];
    check(bao::prepare_density(input, preparation).status() !=
              statistics::DensityStatus::finite,
          "duplicate row identity rejects");
    input = original;
    input.covariance[1] = 2;
    check(bao::prepare_density(input, preparation).status() !=
              statistics::DensityStatus::finite,
          "source asymmetry no repair");
    input = original;
    auto capdp = preparation;
    capdp.maximum_native_bytes = 1;
    auto capped_density = bao::prepare_density(input, capdp);
    check(capped_density.status() != statistics::DensityStatus::finite &&
              capped_density.numerical_status() ==
                  numerics::Status::work_limit &&
              capped_density.source().queries.empty(),
          "source cap before ownedcopy");
    auto cap_eval = dp;
    cap_eval.maximum_native_bytes = 1;
    auto capped_batch = prepared.evaluate(models, cap_eval);
    check(capped_batch.slots.empty() &&
              capped_batch.numerical_status == numerics::Status::work_limit,
          "density output cap before allocate");
    const auto perquery = sizeof(bao::Slot) + sizeof(cosmology::Request) +
                          sizeof(std::size_t) +
                          *cosmology::Expansion::workspace_payload_bound(1);
    auto peak = p;
    peak.maximum_native_bytes = perquery - 1;
    check(bao::evaluate(cq(-1), bao::Ruler(10000),
                        std::span<const bao::Query>(qs.data(), 1), peak)
              .slots.empty(),
          "peak scratch boundary rejects");
    peak.maximum_native_bytes = perquery;
    check(bao::evaluate(cq(-1), bao::Ruler(10000),
                        std::span<const bao::Query>(qs.data(), 1), peak)
                  .slots.size() == 1,
          "exact payload cap boundary admits");
    std::printf("BAO owner checks %d PASS max_budget_fraction %.17g\n", checks,
                worst);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL %s after %d\n", e.what(), checks);
    return 1;
  }
}
