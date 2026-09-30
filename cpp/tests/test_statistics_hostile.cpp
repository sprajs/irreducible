// Independent F03 cases: 2x2 adjugate/rational fixtures, fixed-panel quadrature
// (not production adaptive integration), and a native atanh-series logarithm.
// Named synthetic target budget 1e-10; normalization 1e-9. No generic matrix,
// tail or external math-library qualification follows. No source data required.
#include "../../tests/verification/gaussian_exact_cases.hpp"
#include "irred/statistics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
int checks = 0;
long double max_error = 0;
void check(bool x, const char *why) {
  ++checks;
  if (!x)
    throw std::runtime_error(why);
}
void near(long double actual, long double expected, long double budget,
          const char *why) {
  auto e = std::abs(actual - expected);
  max_error = std::max(max_error, e);
  check(e <= budget, why);
}
constexpr long double ln2pi = 1.8378770664093454835606594728112352797228L;
long double log_reference(long double x) {
  int exponent = 0;
  while (x >= 2) {
    x /= 2;
    ++exponent;
  }
  while (x < 1) {
    x *= 2;
    --exponent;
  }
  auto atanh = [](long double z) {
    long double total = 0, p = z;
    for (int k = 1; k < 200; k += 2) {
      total += p / k;
      p *= z * z;
    }
    return 2 * total;
  };
  return atanh((x - 1) / (x + 1)) + exponent * atanh(1.L / 3);
}
template <class F>
long double simpson(F f, long double a, long double b, unsigned n) {
  long double h = (b - a) / n, sum = f(a) + f(b);
  for (unsigned i = 1; i < n; ++i)
    sum += (i % 2 ? 4 : 2) * f(a + h * i);
  return h * sum / 3;
}
statistics::Metadata metadata() {
  statistics::Metadata m;
  m.ordered_ids = {"A:exposure", "B:exposure"};
  m.measure = "product d(magnitude)";
  m.ordering_provenance = "explicit synthetic row order";
  m.dependence_provenance =
      "unknown source dependence; not an independence assertion";
  return m;
}
} // namespace
int main() {
  using namespace statistics;
  constexpr double budget = 1e-10;
  std::array<double, 4> c{4, 1, 1, 9};
  std::array<double, 2> r{2, -3}, response{1, 1};
  auto m = metadata();
  auto g = prepare_gaussian(c, MatrixKind::covariance, m, 4, budget);
  check(g.status() == DensityStatus::finite, "exact covariance prepared");
  auto d = g.evaluate(r, m.ordered_ids, budget);
  check(d.density.status == DensityStatus::finite, "density finite");
  near(d.quadratic, verifier_gaussian::quadratic.value(), 4e-11,
       "adjugate quadratic");
  near(d.log_determinant, log_reference(35), 4e-11,
       "determinant independent log series");
  near(d.density.log_value, -.5L * (12.L / 5 + log_reference(35) + 2 * ln2pi),
       budget, "normalized Gaussian target");
  std::array<std::string, 2> wrong{m.ordered_ids[1], m.ordered_ids[0]};
  check(g.evaluate(r, wrong, budget).density.status ==
            DensityStatus::incompatible_metadata,
        "measurement order mismatch rejected");
  std::array<double, 4> permuted{9, 1, 1, 4};
  std::array<double, 2> rp{-3, 2};
  auto pm = m;
  pm.ordered_ids.assign(wrong.begin(), wrong.end());
  auto pg = prepare_gaussian(permuted, MatrixKind::covariance, pm, 4, budget);
  near(pg.evaluate(rp, pm.ordered_ids, budget).density.log_value,
       d.density.log_value, budget,
       "joint row covariance permutation invariant");
  std::array<double, 2> zero{0, 0};
  auto density0 = g.evaluate(zero, m.ordered_ids, budget);
  std::array<double, 4> scaled{16, 4, 4, 36};
  auto sg = prepare_gaussian(scaled, MatrixKind::covariance, m, 4, budget);
  near(sg.evaluate(zero, m.ordered_ids, budget).density.log_value -
           density0.density.log_value,
       -log_reference(4), budget,
       "determinant cannot disappear at zero residual");
  std::array<std::size_t, 1> keep{0};
  auto marginal = g.marginal(keep, 4, budget),
       conditional = g.conditional_zero_complement(keep, 4, budget);
  near(marginal.covariance()[0], 4, 4e-11, "marginal variance");
  near(conditional.covariance()[0], 35.L / 9, 4e-11,
       "conditional variance differs");
  std::array<double, 4> w{9. / 35, -1. / 35, -1. / 35, 4. / 35};
  auto precision = prepare_gaussian(w, MatrixKind::precision, m, 4, budget);
  near(precision.marginal(keep, 4, budget).covariance()[0], 4, 4e-11,
       "full precision inversion then marginal selection");
  near(precision.evaluate(r, m.ordered_ids, budget).density.log_value,
       d.density.log_value, budget,
       "precision determinant reciprocal convention");
  check(g.profile_offset(r, response, wrong, budget).status ==
            DensityStatus::incompatible_metadata,
        "profile ordered IDs mismatch rejected");
  check(g.proper_offset(response, wrong, 1. / 3, 2, "wrong-row-prior", true, 4,
                        budget)
                .status() == DensityStatus::incompatible_metadata,
        "proper design vector ordered IDs mismatch rejected");
  auto profile = g.profile_offset(r, response, m.ordered_ids, budget);
  near(profile.coefficient, 7.L / 11, budget, "profile coefficient");
  near(profile.quadratic, 25.L / 11, budget,
       "profile score is not normalized density");
  auto proper = g.proper_offset(response, m.ordered_ids, 1. / 3, 2,
                                "latent-offset-1", true, 4, budget);
  auto p = proper.evaluate(r, m.ordered_ids, budget);
  check(p.density.status == DensityStatus::finite, "proper prior finite");
  near(p.quadratic, 1175.L / 513, 4e-11, "nonzero proper prior mean retained");
  near(p.log_determinant, log_reference(57), 4e-11, "proper determinant");
  near(p.density.log_value,
       -ln2pi - .5L * log_reference(57) - .5L * (1175.L / 513), budget,
       "proper normalized marginal");
  check(
      proper.metadata().dependence_provenance == m.dependence_provenance,
      "declared independent prior does not upgrade unknown source dependence");
  check(proper.proper_offset(response, m.ordered_ids, 1. / 3, 2,
                             "latent-offset-1", true, 4, budget)
                .status() != DensityStatus::finite,
        "latent duplicate refused");
  check(g.proper_offset(response, m.ordered_ids, 1. / 3, 2, "latent-offset-1",
                        false, 4, budget)
                .status() != DensityStatus::finite,
        "latent independence not silently assumed");

  check(proper.priors().size() == 1, "prior record retained");
  const auto &prior = proper.priors()[0];
  near(prior.mean, 1.L / 3, 1e-15, "prior mean public descriptor");
  check(prior.variance == 2 && prior.latent_identity == "latent-offset-1" &&
            prior.independence_assumed,
        "prior variance identity assumption retained");
  check(prior.response == std::vector<double>({1, 1}) &&
            prior.applied_row_ids == m.ordered_ids,
        "original prior response and row IDs retained");
  near(proper.mean_shift()[0], 1.L / 3, 1e-15,
       "prepared mean shift inspectable");
  std::array<double, 2> huge_response{std::numeric_limits<double>::max(), 1};
  check(g.proper_offset(huge_response, m.ordered_ids, 0, 1, "huge", true, 4,
                        budget)
                .status() == DensityStatus::numerical_failure,
        "valid numeric covariance overflow distinguished from invalid input");
  // Integrate conditional likelihood times proper prior directly, using
  // adjugate q(a), not the production covariance update. Posterior
  // N(77/171,70/57):
  // [-16,16] omits <1e-40 mass; fixed panels converge independently.
  auto integrand = [](long double a) {
    auto u = 2 - a, v = -3 - a;
    auto q = (9 * u * u - 2 * u * v + 4 * v * v) / 35;
    auto t = a - 1.L / 3;
    return std::exp(-q / 2 - t * t / 4 - 1.5L * ln2pi -
                    .5L * log_reference(70));
  };
  auto i1 = simpson(integrand, -16, 16, 1024),
       i2 = simpson(integrand, -16, 16, 2048),
       i3 = simpson(integrand, -16, 16, 4096);
  near(i2, i3, 1e-13, "direct prior integral refinement");
  near(i1, i3, 1e-12, "direct prior coarse refinement");
  near(log_reference(i3), p.density.log_value, 1e-9,
       "direct integration normalized proper marginal");
  auto pc = proper.conditional_zero_complement(keep, 4, budget);
  std::array<double, 1> ra{2};
  near(pc.covariance()[0], 57.L / 11, 4e-11, "proper conditional covariance");
  near(pc.evaluate(ra, pc.metadata().ordered_ids, budget).quadratic,
       275.L / 513, 4e-11,
       "conditional complement centered zero with retained prior mean");

  check(pc.priors()[0].applied_row_ids == m.ordered_ids &&
            pc.priors()[0].response == prior.response,
        "selection preserves original prior ancestry");
  check(pc.selection_history().size() == 1 &&
            pc.selection_history()[0].kept_row_ids ==
                std::vector<std::string>{m.ordered_ids[0]} &&
            pc.selection_history()[0].complement_row_ids ==
                std::vector<std::string>{m.ordered_ids[1]},
        "conditional selection identity retained");
  check(selected_standard_normal_positive(0).status ==
            DensityStatus::outside_support,
        "strict support boundary");
  check(selected_standard_normal_positive(-1).status ==
            DensityStatus::outside_support,
        "negative selected support");
  check(selected_standard_normal_positive(INFINITY).status ==
            DensityStatus::invalid_input,
        "nonfinite is failure not support");
  near(selected_standard_normal_positive(1).log_value -
           normal_log_density(1, 0, 1).log_value,
       log_reference(2), budget, "selection denominator");
  // A singleton endpoint is measure zero; use right-hand limiting density in
  // quadrature, while the public support test at exactly zero remains strict.
  auto half = [](long double x) {
    if (x == 0)
      return std::exp(-ln2pi / 2) * 2;
    return std::exp(static_cast<long double>(
        selected_standard_normal_positive(static_cast<double>(x)).log_value));
  };
  near(simpson(half, 0, 8, 2048), 1, 1e-9,
       "selected density normalizes; tail <1.3e-15");
  near(normal_log_density(2, 0, 2).log_value -
           normal_log_density(1, 0, 1).log_value,
       -log_reference(2), budget, "normal density measure unit scaling");
  check(normal_log_density(0, 0, -1).status == DensityStatus::invalid_input,
        "invalid normal variance");
  check(normal_log_density(std::numeric_limits<double>::max(), 0,
                           std::numeric_limits<double>::min())
                .status == DensityStatus::numerical_failure,
        "finite normal overflow not outside support");
  check(poisson_log_mass(0, 0).status == DensityStatus::finite &&
            poisson_log_mass(0, 0).log_value == 0,
        "zero rate zero count mass one");
  check(poisson_log_mass(1, 0).status == DensityStatus::outside_support,
        "zero rate positive count exact support");
  check(poisson_log_mass(100, 1).status == DensityStatus::unsupported_domain,
        "gamma count domain respected");
  near(poisson_log_mass(2, 3).log_value,
       2 * log_reference(3) - 3 - log_reference(2), budget,
       "Poisson independent factorial");
  auto bad = m;
  bad.ordered_ids[1] = bad.ordered_ids[0];
  check(prepare_gaussian(c, MatrixKind::covariance, bad, 4, budget).status() ==
            DensityStatus::incompatible_metadata,
        "duplicate measurements rejected");
  std::array<double, 4> asym{4, 1.00000003, 1, 9};
  check(prepare_gaussian(asym, MatrixKind::covariance, m, 4, budget).status() ==
            DensityStatus::numerical_failure,
        "raw source asymmetry never repaired");
  std::array<double, 4> singular{1, 1, 1, 1};
  check(prepare_gaussian(singular, MatrixKind::covariance, m, 4, budget)
                .status() == DensityStatus::numerical_failure,
        "singular no jitter");
  std::array<double, 2> nan{NAN, 0};
  check(g.evaluate(nan, m.ordered_ids, budget).density.status ==
            DensityStatus::numerical_failure,
        "nonfinite residual failure");
  check(g.evaluate(r, m.ordered_ids, 0).density.status ==
            DensityStatus::numerical_failure,
        "solve consumer budget required");

  // Rational Householder rotation of diag(1,4,16), determinant64.
  // Adjugate inverse is exact rational below; separate matrix construction
  // challenges derived inverse symmetry and selection at normalized-target
  // level.
  std::array<double, 9> c3{41. / 9,  52. / 9, 4. / 9,   52. / 9, 116. / 9,
                           -16. / 9, 4. / 9,  -16. / 9, 32. / 9};
  std::array<double, 9> w3{2. / 3, -1. / 3, -1. / 4, -1. / 3, 1. / 4,
                           1. / 6, -1. / 4, 1. / 6,  19. / 48};
  std::array<double, 3> r3{2, -1, 3};
  auto m3 = m;
  m3.ordered_ids.push_back("C:exposure");
  auto g3 = prepare_gaussian(c3, MatrixKind::covariance, m3, 9, budget),
       pw3 = prepare_gaussian(w3, MatrixKind::precision, m3, 9, budget);
  check(g3.status() == DensityStatus::finite &&
            pw3.status() == DensityStatus::finite,
        "rotated3 covariance and precision prepared");
  auto v3 = g3.evaluate(r3, m3.ordered_ids, budget),
       vp3 = pw3.evaluate(r3, m3.ordered_ids, budget);
  near(v3.quadratic, 61.L / 16, 4e-11, "rational rotated3 quadratic");
  near(v3.log_determinant, log_reference(64), 4e-11, "rotated3 determinant");
  near(v3.density.log_value, -.5L * (61.L / 16 + log_reference(64) + 3 * ln2pi),
       budget, "rotated3 normalized density");
  near(vp3.density.log_value, v3.density.log_value, budget,
       "precision derived operator normalized parity3");
  check(v3.backward_residual < 1e-12 && vp3.backward_residual < 1e-12,
        "assembled3 solve residual budgets");
  std::array<std::size_t, 2> k3{0, 2};
  std::array<double, 2> rk3{2, 3};
  auto mg3 = pw3.marginal(k3, 9, budget),
       cg3 = pw3.conditional_zero_complement(k3, 9, budget);
  near(mg3.evaluate(rk3, mg3.metadata().ordered_ids, budget).density.log_value,
       -.5L * (449.L / 144 + log_reference(16) + 2 * ln2pi), budget,
       "precision3 marginal normalized target");
  near(cg3.evaluate(rk3, cg3.metadata().ordered_ids, budget).density.log_value,
       -.5L * (155.L / 48 + log_reference(144.L / 29) + 2 * ln2pi), budget,
       "precision3 centered conditional normalized target differs");
  std::printf("{\"suite\":\"independent_statistics\",\"checks\":%d,\"maximum_"
              "error\":%.17Lg,\"passed\":true}\n",
              checks, max_error);
}
