#include "irred/correlated_calibration.hpp"
#include "irred/gaussian_posterior.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred::statistics;
namespace {
unsigned checks = 0;
void need(bool yes, const char *why) {
  ++checks;
  if (!yes)
    throw std::runtime_error(why);
}
void near(long double x, long double y) {
  need(std::abs(x - y) <= 2e-12L * (1 + std::abs(y)),
       "frozen independent allocation");
}
Metadata metadata() {
  Metadata m;
  m.ordered_ids = {"r0", "r1"};
  m.measure = "dr0 dr1";
  m.table_identity = "synthetic fixed Gaussian";
  m.ordering_provenance = "explicit";
  return m;
}
Gaussian source() {
  return prepare_gaussian(std::array<double, 4>{2, .25, .25, 1.5},
                          MatrixKind::covariance, metadata(), 100, 1e-10,
                          irred::numerics::Arithmetic::longdouble_cpu_v1);
}
ParameterPrior prior() {
  return {{"zero", "colour"},
          {"mag", "mag"},
          {"zero"},
          {.5, -.25},
          {1, .375, .375, .5},
          "synthetic proper prior",
          "explicit fixed X",
          "mag",
          "d(zero) d(colour)",
          "explicit noise/prior independence",
          true};
}
} // namespace
int main() {
  try {
    const std::array<double, 4> x{1, .5, -.25, 1};
    const std::array<double, 2> r{1.25, -.5};
    auto ids = metadata().ordered_ids;
    auto p = prior();
    auto g = source();
    const auto bound = GaussianPosterior::preparation_payload_bound(g, p);
    need(bound.has_value(), "payload exists");
    auto owner = GaussianPosterior::prepare(std::move(g), x, ids, p);
    need(owner.status() == DensityStatus::finite &&
             g.status() == DensityStatus::invalid_input,
         "owner and consumed source");
    const auto mean = owner.condition(r, ids);
    need(mean.status == DensityStatus::finite && mean.value.size() == 2,
         "mean admitted");
    // Exact rational cofactor oracle, independently derived before
    // implementation.
    near(mean.value[0], 20604.L / 25511);
    near(mean.value[1], -7125.L / 51022);
    const long double v[]{15160.L / 25511, 4466.L / 25511, 4466.L / 25511,
                          8592.L / 25511};
    for (unsigned j = 0; j < 4; ++j)
      near(owner.covariance()[j], v[j]);
    const std::array<double, 2> beta{.7, -.2};
    auto density = owner.log_density(r, ids, beta, p.ordered_parameter_ids);
    need(density.density.status == DensityStatus::finite,
         "normalized parameter density admitted");
    const long double a = beta[0] - 20604.L / 25511,
                      b = beta[1] + 7125.L / 51022;
    const long double q = (2148.L / 1081) * a * a -
                          2 * (2233.L / 2162) * a * b + (3790.L / 1081) * b * b;
    near(density.quadratic, q);
    near(density.log_determinant, std::log(v[0] * v[3] - v[1] * v[2]));
    near(density.density.log_value,
         -.5L * (q + std::log(v[0] * v[3] - v[1] * v[2]) +
                 2 * std::log(2 * std::numbers::pi_v<long double>)));
    // Independently composed Bayes identity: shared factor kernels are
    // ancestry.
    auto prior_g = prepare_gaussian(p.covariance, MatrixKind::covariance,
                                    owner.source().metadata(), 100, 1e-10,
                                    owner.source().arithmetic());
    std::array<double, 2> prior_r{beta[0] - .5, beta[1] + .25},
        noise_r{r[0] - beta[0] - .5 * beta[1], r[1] + .25 * beta[0] - beta[1]};
    CalibrationPrior cp{p.ordered_parameter_ids,
                        p.parameter_units,
                        p.mean,
                        p.covariance,
                        p.prior_identity,
                        p.design_identity,
                        p.residual_unit,
                        "synthetic",
                        p.dependence_identity,
                        metadata().measure,
                        true};
    auto sg = source();
    auto predictive = CorrelatedCalibration::prepare(std::move(sg), x, ids, cp);
    const auto evidence = predictive.evaluate(r, ids);
    near(density.density.log_value,
         owner.source().evaluate(noise_r, ids, 1e-10).density.log_value +
             prior_g.evaluate(prior_r, ids, 1e-10).density.log_value -
             evidence.density.log_value);
    // Marginally unobserved directions remain proper and prior-correlated.
    auto zero_source = source();
    auto zero = GaussianPosterior::prepare(
        std::move(zero_source), std::array<double, 4>{0, 0, 0, 0}, ids, p);
    need(zero.status() == DensityStatus::finite,
         "null design admitted by proper prior");
    auto zm = zero.condition(r, ids);
    near(zm.value[0], p.mean[0]);
    near(zm.value[1], p.mean[1]);
    for (unsigned j = 0; j < 4; ++j)
      near(zero.covariance()[j], p.covariance[j]);
    auto rank_source = source();
    auto rank = GaussianPosterior::prepare(
        std::move(rank_source), std::array<double, 4>{1, 1, 2, 2}, ids, p);
    need(rank.status() == DensityStatus::finite &&
             rank.condition(r, ids).status == DensityStatus::finite,
         "rank-deficient likelihood plus proper prior");
    // Parameter-coordinate change has its density Jacobian, distinct from
    // units.
    auto scaled = p;
    scaled.mean[0] *= 10;
    scaled.covariance[0] *= 100;
    scaled.covariance[1] *= 10;
    scaled.covariance[2] *= 10;
    scaled.parameter_units[0] = "0.1mag";
    scaled.parameter_measure = "d(scaled_zero) d(colour)";
    auto scaled_source = source();
    auto so = GaussianPosterior::prepare(
        std::move(scaled_source), std::array<double, 4>{.1, .5, -.025, 1}, ids,
        scaled);
    need(so.status() == DensityStatus::finite,
         "unit-transformed posterior admitted");
    auto sm = so.condition(r, ids);
    near(sm.value[0], 10 * mean.value[0]);
    near(sm.value[1], mean.value[1]);
    near(so.log_density(r, ids, std::array<double, 2>{7, -.2},
                        scaled.ordered_parameter_ids)
             .density.log_value,
         density.density.log_value - std::log(10.L));
    for (unsigned mode = 0; mode < 7; ++mode) {
      auto bad = p;
      auto original = source();
      PosteriorPolicy policy;
      if (mode == 0)
        bad.noise_independence_declared = false;
      if (mode == 1)
        bad.covariance = {1, 1, 1, 1};
      if (mode == 2)
        bad.parameter_measure.clear();
      if (mode == 3)
        bad.ordered_parameter_ids[1] = bad.ordered_parameter_ids[0];
      if (mode == 4)
        bad.shared_nuisance_ids = {"missing"};
      if (mode == 5)
        policy.maximum_payload_bytes = 0;
      if (mode == 6)
        policy.maximum_work_units = 0;
      auto failed =
          GaussianPosterior::prepare(std::move(original), x, ids, bad, policy);
      need(failed.status() != DensityStatus::finite &&
               original.status() == DensityStatus::finite,
           "failed preparation retains source");
    }
    need(owner.condition(r, std::array<std::string, 2>{"r1", "r0"}).status ==
             DensityStatus::incompatible_metadata,
         "ordered rows exact");
    need(owner.log_density(r, ids, beta,
                           std::array<std::string, 2>{"colour", "zero"})
                 .density.status == DensityStatus::incompatible_metadata,
         "ordered parameters exact");
    PosteriorPolicy impossible;
    impossible.maximum_forward_sensitivity = 1e-30;
    need(owner.condition(r, ids, impossible).status != DensityStatus::finite,
         "unattainable mean budget refused");
    PosteriorPolicy tiny;
    tiny.maximum_payload_bytes = 0;
    need(owner.log_density(r, ids, beta, p.ordered_parameter_ids, tiny)
                 .density.status != DensityStatus::finite,
         "evaluation quota admitted before vectors");
    auto moved = std::move(owner);
    need(owner.status() == DensityStatus::invalid_input &&
             owner.condition(r, ids).value.empty(),
         "move invalidates source");
    auto *alias = &moved;
    moved = std::move(*alias);
    near(moved.condition(r, ids).value[0], mean.value[0]);
    GaussianPosterior assigned;
    assigned = std::move(moved);
    need(moved.status() == DensityStatus::invalid_input,
         "move assignment invalidates source");
    near(assigned.condition(r, ids).value[1], mean.value[1]);
    int previous = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    auto unsupported = assigned.condition(r, ids);
    std::fesetround(previous);
    need(unsupported.status != DensityStatus::finite,
         "unsupported rounding refused");
    std::cout << "PASS " << checks << " Gaussian posterior controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
