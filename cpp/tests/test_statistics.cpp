// Synthetic F03 acceptance only. Predeclared log-density budget 1e-10;
// q/logdet/normalization component budgets 4e-11/4e-11/2e-11.
// Independent rational references: C=[4,1;1,9], det35, inverse=[9,-1;-1,4]/35.
// Direct proper-prior integration uses scalar product densities, not covariance
// assembly. Normalization budget1e-9, analytic truncation tails <1e-10.
// Logarithms in references share system libm: not independent libm
// qualification.
#include "irred/statistics.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred;
namespace {
int checks = 0;
double maximum = 0;
void check(bool v, const char *n) {
  ++checks;
  if (!v)
    throw std::runtime_error(n);
}
void near(double v, long double ref, double b, const char *n) {
  auto e = std::abs(static_cast<long double>(v) - ref);
  maximum = std::max(maximum, static_cast<double>(e));
  check(e <= b, n);
}
statistics::Metadata meta() {
  statistics::Metadata m;
  m.ordered_ids = {"a", "b"};
  m.measure = "dmag_a dmag_b";
  m.ordering_provenance = "explicit synthetic ordered IDs";
  return m;
}
double prior_integrand(double a, const void *) {
  const long double pi = std::numbers::pi_v<long double>;
  return static_cast<double>(
      std::exp(-(2 - a) * (2 - a) / 6 - (-3 - a) * (-3 - a) / 16 - a * a / 2) /
      (std::pow(2 * pi, 1.5L) * std::sqrt(24.L)));
}
// The boundary has measure zero; quadrature uses its right-hand limiting
// density to avoid an artificial endpoint jump while support tests retain
// strict x>0.
double halfnormal(double x, const void *) {
  if (x == 0)
    return std::sqrt(2 / std::numbers::pi);
  auto r = statistics::selected_standard_normal_positive(x);
  return r.status == statistics::DensityStatus::finite ? std::exp(r.log_value)
                                                       : 0;
}
double jacobian(double u, const void *) { return std::exp(u - std::exp(u)); }
double normal(double x, const void *) {
  return std::exp(statistics::normal_log_density(x, 0, 1).log_value);
}
} // namespace
int main() {
  try {
    using namespace statistics;
    constexpr double budget = 1e-10;
    const auto pi = std::numbers::pi_v<long double>;
    near(normal_log_density(0, 0, 1).log_value, -std::log(2 * pi) / 2, 2e-11,
         "normal normalized center");
    near(normal_log_density(2, 0, 2).log_value,
         -.5L - std::log(2.L) - std::log(2 * pi) / 2, budget,
         "normal scale determinant");
    check(normal_log_density(0, 0, 0).status == DensityStatus::invalid_input,
          "zero normal scale");
    check(normal_log_density(INFINITY, 0, 1).status ==
              DensityStatus::invalid_input,
          "nonfinite normal");
    check(normal_log_density(std::numeric_limits<double>::max(),
                             -std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::min())
                  .status == DensityStatus::numerical_failure,
          "finite normal overflow not support");
    check(poisson_log_mass(0, 0).status == DensityStatus::finite &&
              poisson_log_mass(0, 0).log_value == 0,
          "zero count zero rate");
    check(poisson_log_mass(1, 0).status == DensityStatus::outside_support,
          "zero rate positive count");
    check(poisson_log_mass(0, -1).status == DensityStatus::invalid_input,
          "negative rate");
    check(poisson_log_mass(100, 2).status == DensityStatus::unsupported_domain,
          "count beyond gamma scope");
    near(poisson_log_mass(2, 3).log_value,
         2 * std::log(3.L) - 3 - std::log(2.L), budget, "factorial mass");
    long double sum = 0;
    for (unsigned k = 0; k < 40; ++k)
      sum +=
          std::exp(static_cast<long double>(poisson_log_mass(k, 3).log_value));
    near(static_cast<double>(sum), 1, 1e-12,
         "Poisson normalization tail geometric <1e-29");
    std::array<double, 4> c{4, 1, 1, 9};
    std::array<double, 2> r{2, -3};
    auto g = prepare_gaussian(c, MatrixKind::covariance, meta(), 4, budget);
    check(g.status() == DensityStatus::finite, "Gaussian preparation");
    auto v = g.evaluate(r, g.metadata().ordered_ids, budget);
    check(v.density.status == DensityStatus::finite, "Gaussian evaluated");
    near(v.quadratic, 12.L / 5, 4e-11, "rational quadratic");
    near(v.log_determinant, std::log(35.L), 4e-11, "exact determinant");
    near(v.normalization, 2 * std::log(2 * pi), 2e-11, "measure normalization");
    near(v.density.log_value,
         -.5L * (12.L / 5 + std::log(35.L) + 2 * std::log(2 * pi)), budget,
         "assembled density");
    auto wrong = meta().ordered_ids;
    std::swap(wrong[0], wrong[1]);
    check(g.evaluate(r, wrong, budget).density.status ==
              DensityStatus::incompatible_metadata,
          "order mismatch");
    std::array<std::size_t, 1> keep{0};
    near(g.marginal(keep, 4, budget).covariance()[0], 4, 4e-11,
         "marginal covariance");
    near(g.conditional_zero_complement(keep, 4, budget).covariance()[0],
         35.L / 9, 4e-11, "conditional covariance");
    std::array<double, 4> w{9. / 35, -1. / 35, -1. / 35, 4. / 35};
    auto precision =
        prepare_gaussian(w, MatrixKind::precision, meta(), 4, budget);
    near(precision.marginal(keep, 4, budget).covariance()[0], 4, 4e-11,
         "precision marginal uses full inverse action");
    near(precision.evaluate(r, precision.metadata().ordered_ids, budget)
             .density.log_value,
         v.density.log_value, budget, "precision normalized parity");
    std::array<double, 2> x{1, 1};
    check(g.profile_offset(r, x, wrong, budget).status ==
              DensityStatus::incompatible_metadata,
          "wrong profile row order rejected");
    check(g.proper_offset(x, wrong, 0, 1, "wrong-design", true, 4, budget)
                  .status() == DensityStatus::incompatible_metadata,
          "wrong prior response row order rejected");
    auto profile = g.profile_offset(r, x, meta().ordered_ids, budget);
    near(profile.coefficient, 7.L / 11, budget, "offset profile coefficient");
    near(profile.quadratic, 25.L / 11, budget,
         "offset profile rational projected q");
    std::array<double, 4> d{3, 0, 0, 8};
    auto noise = prepare_gaussian(d, MatrixKind::covariance, meta(), 4, budget);
    auto marginal = noise.proper_offset(x, meta().ordered_ids, 0, 1, "offset-A",
                                        true, 4, budget);
    near(marginal.evaluate(r, marginal.metadata().ordered_ids, budget)
             .density.log_value,
         v.density.log_value, budget, "proper nuisance covariance");
    check(marginal.proper_offset(x, meta().ordered_ids, 0, 1, "offset-A", true,
                                 4, budget)
                  .status() != DensityStatus::finite,
          "same latent identity cannot be counted twice");
    check(noise.proper_offset(x, meta().ordered_ids, 0, 1, "offset-A", false, 4,
                              budget)
                  .status() != DensityStatus::finite,
          "unknown latent independence rejected");
    for (double tol : {1e-8, 1e-10, 1e-12}) {
      auto integral = numerics::integrate(prior_integrand, nullptr, -10, 10,
                                          {tol, 0, 100000, 30});
      check(integral.status == numerics::Status::ok,
            "prior integral converges");
      if (tol == 1e-12)
        near(std::log(integral.value), v.density.log_value, 1e-9,
             "independent proper integration posterior mean1/5 variance24/35 "
             "tails negligible");
    }
    auto nonzero = g.proper_offset(x, meta().ordered_ids, 1. / 3, 2,
                                   "nonzero-offset", true, 4, budget);
    auto nz = nonzero.evaluate(r, nonzero.metadata().ordered_ids, budget);
    near(nz.log_determinant, std::log(57.L), 4e-11,
         "nonzero proper determinant");
    near(nz.quadratic, 1175.L / 513, 4e-11, "nonzero prior mean included");
    check(nonzero.priors().size() == 1 && nonzero.priors()[0].mean == 1. / 3 &&
              nonzero.priors()[0].variance == 2 &&
              nonzero.priors()[0].response == std::vector<double>({1, 1}) &&
              nonzero.priors()[0].applied_row_ids == meta().ordered_ids &&
              nonzero.priors()[0].independence_assumed,
          "owned reconstructible original prior");
    check(nonzero.metadata().dependence_provenance.empty(),
          "assumed prior independence does not invent source lineage");
    auto reduced = nonzero.conditional_zero_complement(keep, 4, budget);
    check(reduced.priors()[0].applied_row_ids.size() == 2 &&
              reduced.selection_history().size() == 1 &&
              reduced.selection_history()[0].complement_row_ids ==
                  std::vector<std::string>({"b"}),
          "original prior not mislabeled conditional posterior");
    near(reduced.mean_shift()[0], 1.L / 3, budget,
         "centered conditional mean retained");
    check(g.proper_offset(std::array<double, 2>{2, 2}, meta().ordered_ids, 0,
                          std::numeric_limits<double>::max(), "overflow", true,
                          4, budget)
                  .status() == DensityStatus::numerical_failure,
          "proper arithmetic overflow is numerical failure");
    std::array<double, 2> shared{1, 2};
    auto joint = noise.proper_offset(shared, meta().ordered_ids, 0, 1,
                                     "shared-calibration", true, 4, budget);
    auto j = joint.evaluate(r, joint.metadata().ordered_ids, budget);
    near(j.quadratic, 27.L / 11, 4e-11, "shared calibration q");
    near(j.log_determinant, std::log(44.L), 4e-11, "shared calibration det");
    check(joint.covariance()[1] == 2, "shared crossblock retained");
    std::array<double, 4> separate{4, 0, 0, 12};
    auto falsely_independent =
        prepare_gaussian(separate, MatrixKind::covariance, meta(), 4, budget);
    check(std::abs(falsely_independent.evaluate(r, meta().ordered_ids, budget)
                       .density.log_value -
                   j.density.log_value) > 1e-3,
          "discarded crossblock changes density");
    check(selected_standard_normal_positive(0).status ==
              DensityStatus::outside_support,
          "strict selected boundary");
    near(selected_standard_normal_positive(1).log_value -
             normal_log_density(1, 0, 1).log_value,
         std::log(2.L), budget, "selection divisor");
    auto half =
        numerics::integrate(halfnormal, nullptr, 0, 8, {1e-12, 0, 100000, 30});
    check(half.status == numerics::Status::ok, "selected integral converges");
    near(half.value, 1, 1e-9, "selected normalization tail<1.3e-15");
    auto norm =
        numerics::integrate(normal, nullptr, -8, 8, {1e-12, 0, 100000, 30});
    near(norm.value, 1, 1e-9, "normal normalization");
    auto jac =
        numerics::integrate(jacobian, nullptr, -26, 4, {1e-12, 0, 100000, 30});
    check(jac.status == numerics::Status::ok, "Jacobian integral converges");
    near(jac.value, 1, 1e-9,
         "data Jacobian normalization omitted mass<5.2e-12");
    std::array<double, 4> singular{1, 1, 1, 1}, bad{1, 2, 2, 1};
    check(prepare_gaussian(singular, MatrixKind::covariance, meta(), 4, budget)
                  .status() == DensityStatus::numerical_failure,
          "singular no jitter");
    check(prepare_gaussian(bad, MatrixKind::covariance, meta(), 4, budget)
                  .status() == DensityStatus::numerical_failure,
          "indefinite rejects");
    check(g.evaluate(r, g.metadata().ordered_ids, 1e-30).density.status ==
              DensityStatus::numerical_failure,
          "conditioning policy honored");
    observations::Input input{};
    input.profile = observations::Profile::gaussian_fixture_v1;
    input.role = observations::Role::synthetic_control;
    input.unit = observations::Unit::magnitude;
    input.calibration = observations::Calibration::unknown;
    input.uncertainty = observations::Uncertainty::precision;
    input.uncertainty_unit =
        observations::UncertaintyUnit::inverse_magnitude_squared;
    input.table_sha256 = std::string(64, 'a');
    input.uncertainty_sha256 = std::string(64, 'b');
    input.ordering_provenance = "explicit synthetic source axes";
    input.measurement_ids = {"a", "b"};
    input.event_ids = {"same-event", "same-event"};
    input.uncertainty_axis_ids = input.measurement_ids;
    input.values = {2, -3};
    input.missing = {0, 0};
    input.quality = {0, 0};
    input.source_selection = {1, 0};
    input.uncertainty_matrix.assign(w.begin(), w.end());
    auto observation = observations::prepare(std::move(input), {2, 4});
    auto consumer = prepare_observations(
        observation, observations::Selection::all, 4, budget);
    check(consumer.status() == DensityStatus::finite,
          "typed observation consumer prepares");
    near(consumer.covariance()[0], 4, 4e-11,
         "selected precision is marginal not conditional");
    check(consumer.metadata().ordered_ids == std::vector<std::string>({"a"}) &&
              consumer.metadata().table_identity == std::string(64, 'a') &&
              consumer.metadata().uncertainty_identity == std::string(64, 'b'),
          "source identities and selected ordering retained");
    check(consumer.metadata().calibration_provenance.empty() &&
              consumer.metadata().dependence_provenance.empty() &&
              !consumer.metadata().source_semantics.empty(),
          "unknown lineage remains unknown in consumer");
    std::printf("{\"suite\":\"statistics_contract\",\"checks\":%d,\"max_"
                "absolute_error\":%.17g,\"passed\":true}\n",
                checks, maximum);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
