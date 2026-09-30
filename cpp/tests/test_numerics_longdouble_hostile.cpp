// Independent named higher-precision checks: dyadic exact solutions, rational
// 2x2 adjugate, canonical binary64 3x3 cofactor inverse, and a rotated dyadic
// spectrum. No production Cholesky/LDLT ancestry in expected solutions.
// Absolute/relative solution budget1e-10 for these named small fixtures only.
#include "irred/numerics.hpp"
#include "irred/statistics.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred::numerics;
namespace {
int checks = 0;
void check(bool ok, const char *why) {
  ++checks;
  if (!ok)
    throw std::runtime_error(why);
}
void near(long double actual, long double expected, long double abs,
          long double rel, const char *why) {
  check(std::abs(actual - expected) <= abs + rel * std::abs(expected), why);
}
struct RoundingRestore {
  std::fenv_t environment{};
  RoundingRestore() { std::fegetenv(&environment); }
  ~RoundingRestore() { std::fesetenv(&environment); }
};
std::array<long double, 3> cofactor_solution(const std::array<double, 9> &c,
                                             const std::array<double, 3> &b) {
  const long double a = c[0], d = c[4], f = c[8], x = c[1], y = c[2], z = c[5];
  const auto det =
      a * (d * f - z * z) - x * (x * f - y * z) + y * (x * z - y * d);
  return {((d * f - z * z) * b[0] + (y * z - x * f) * b[1] +
           (x * z - y * d) * b[2]) /
              det,
          ((y * z - x * f) * b[0] + (a * f - y * y) * b[1] +
           (x * y - a * z) * b[2]) /
              det,
          ((x * z - y * d) * b[0] + (x * y - a * z) * b[1] +
           (a * d - x * x) * b[2]) /
              det};
}
} // namespace
int main() {
  RoundingRestore restore;
  check(std::fesetround(FE_TONEAREST) == 0, "rounding control available");
  check(std::numeric_limits<long double>::radix == 2 &&
            std::numeric_limits<long double>::digits >= 64,
        "qualified host longdouble arithmetic");
  constexpr auto wide = Arithmetic::longdouble_cpu_v1;
  std::array<double, 4> diagonal{4, 0, 0, 16};
  std::array<double, 2> rhs{1, 4};
  auto implicit = cholesky(diagonal, 2, 4),
       explicit_legacy =
           cholesky(diagonal, 2, 4, Arithmetic::binary64_legacy_v1),
       ld = cholesky(diagonal, 2, 4, wide);
  check(ld.status() == Status::ok && ld.arithmetic() == wide,
        "explicit prepared arithmetic identity");
  auto baseline = solve(implicit, rhs, 1e-10),
       same = solve(explicit_legacy, rhs, 1e-10), exact = solve(ld, rhs, 1e-10);
  check(baseline.value == same.value &&
            baseline.backward_residual == same.backward_residual &&
            baseline.estimated_forward_sensitivity ==
                same.estimated_forward_sensitivity,
        "default legacy remains unchanged");
  check(exact.status == Status::ok &&
            exact.value == std::vector<double>({.25, .25}) &&
            exact.output_rounding_error_inf == 0,
        "dyadic exact finalcast");
  std::array<double, 4> c{4, 1, 1, 9};
  std::array<double, 2> ones{1, 1}, zero{0, 0};
  auto rational = cholesky(c, 2, 4, wide);
  auto x = solve(rational, ones, 1e-10);
  check(x.status == Status::ok, "rational wider solve");
  near(x.value[0], 8.L / 35, 1e-16, 0, "adjugate rational component0");
  near(x.value[1], 3.L / 35, 1e-16, 0, "adjugate rational component1");
  check(x.output_rounding_error_inf > 0 && x.output_rounding_error_relative > 0,
        "binary64 cast error retained separately");
  const auto residual0 = 4.L * x.value[0] + x.value[1] - 1,
             residual1 = x.value[0] + 9.L * x.value[1] - 1;
  const auto expected_backward =
      std::max(std::abs(residual0), std::abs(residual1)) /
      (10.L * std::max(std::abs(x.value[0]), std::abs(x.value[1])) + 1);
  near(x.post_cast_backward_residual, expected_backward, 1e-30, 1e-12,
       "diagnostic refers to returned vector");
  check(x.backward_residual == x.post_cast_backward_residual,
        "public residual names agree");
  auto zeros = solve(rational, zero, 1e-10);
  check(zeros.status == Status::ok &&
            zeros.value == std::vector<double>({0, 0}) &&
            zeros.output_rounding_error_inf == 0,
        "zero solve/cast policy");
  std::array<double, 9> rotated{41. / 9,  52. / 9, 4. / 9,   52. / 9, 116. / 9,
                                -16. / 9, 4. / 9,  -16. / 9, 32. / 9};
  std::array<double, 3> b{2, -1, 3};
  auto f3 = cholesky(rotated, 3, 9, wide);
  auto s3 = solve(f3, b, 1e-10);
  auto reference = cofactor_solution(rotated, b);
  check(s3.status == Status::ok, "canonical binary64 rotated3 solve");
  for (std::size_t i = 0; i < 3; ++i)
    near(s3.value[i], reference[i], 1e-14, 1e-10,
         "independent cofactor inverse");
  const double delta = std::ldexp(1., -20);
  std::array<double, 4> ill{(1 + delta) / 2, (1 - delta) / 2, (1 - delta) / 2,
                            (1 + delta) / 2};
  std::array<double, 2> anti{1, -1};
  auto ill_old = cholesky(ill, 2, 4), ill_wide = cholesky(ill, 2, 4, wide);
  check(solve(ill_old, anti, 1e-10).status ==
            Status::conditioning_budget_exceeded,
        "legacy conservative guard still rejects");
  auto ill_solution = solve(ill_wide, anti, 1e-10);
  check(ill_solution.status == Status::ok,
        "named wider arithmetic meets unchanged requested guard");
  near(ill_solution.value[0], 1.L / delta, 0, 1e-10,
       "rotated dyadic eigensystem0");
  near(ill_solution.value[1], -1.L / delta, 0, 1e-10,
       "rotated dyadic eigensystem1");
  std::array<double, 4> near_singular{1, 0, 0, 1e-16};
  std::array<double, 2> all_one{1, 1};
  check(solve(cholesky(near_singular, 2, 4, wide), all_one, 1e-10).status ==
            Status::conditioning_budget_exceeded,
        "higher precision does not waive ill conditioning");
  std::array<double, 1> large_diagonal{1e6},
      tiny_rhs{std::numeric_limits<double>::min()}, tiny_diagonal{1e-200},
      large_rhs{1e200};
  check(solve(cholesky(large_diagonal, 1, 1, wide), tiny_rhs, 1e-10).status ==
            Status::outside_domain,
        "nonzero returned subnormal fails explicit domain");
  check(solve(cholesky(tiny_diagonal, 1, 1, wide), large_rhs, 1e-10).status ==
            Status::overflow,
        "wide work cannot hide binary64 output overflow");
  check(cholesky(c, 2, 4, static_cast<Arithmetic>(99)).status() ==
            Status::invalid_input,
        "unknown arithmetic rejected");
  check(std::fesetround(FE_UPWARD) == 0, "alternate rounding available");
  check(cholesky(c, 2, 4, wide).status() == Status::outside_domain,
        "prepare rejects unsupported rounding");
  check(solve(rational, ones, 1e-10).status == Status::outside_domain,
        "solve rechecks changed environment");
  std::fesetround(FE_TONEAREST);
  std::feclearexcept(FE_ALL_EXCEPT);
  std::feraiseexcept(FE_INVALID);
  const auto sticky = std::fetestexcept(FE_ALL_EXCEPT);
  auto sticky_factor = cholesky(c, 2, 4, wide);
  check(sticky_factor.status() == Status::ok &&
            std::fetestexcept(FE_ALL_EXCEPT) == sticky,
        "prepare does not mistake or discard caller sticky flags");
  auto sticky_solution = solve(sticky_factor, ones, 1e-10);
  check(sticky_solution.status == Status::ok &&
            std::fetestexcept(FE_ALL_EXCEPT) == sticky,
        "solve restores caller exception flags");
  std::feclearexcept(FE_ALL_EXCEPT);
  irred::statistics::Metadata metadata;
  metadata.ordered_ids = {"synthetic:A", "synthetic:B"};
  metadata.measure = "product d(magnitude)";
  metadata.ordering_provenance = "explicit canonical fixture axes";
  auto gaussian = irred::statistics::prepare_gaussian(
      c, irred::statistics::MatrixKind::covariance, metadata, 4, 1e-10, wide);
  check(gaussian.status() == irred::statistics::DensityStatus::finite &&
            gaussian.arithmetic() == wide &&
            gaussian.metadata().arithmetic_id == "F02/longdouble-cpu/v1",
        "statistics records actual arithmetic");
  std::array<std::size_t, 1> kept{0};
  auto marginal = gaussian.marginal(kept, 1, 1e-10);
  auto conditional = gaussian.conditional_zero_complement(kept, 1, 1e-10);
  check(marginal.arithmetic() == wide && conditional.arithmetic() == wide,
        "subset policies inherit prepared arithmetic");
  auto proper = gaussian.proper_offset(ones, metadata.ordered_ids, 1. / 3, 2,
                                       "synthetic-latent", true, 4, 1e-10);
  check(proper.status() == irred::statistics::DensityStatus::finite &&
            proper.arithmetic() == wide,
        "proper covariance update inherits arithmetic");
  auto profile = std::move(gaussian).prepare_offset_profile(
      ones, metadata.ordered_ids, 1e-10);
  check(profile.status() == irred::statistics::DensityStatus::finite &&
            profile.arithmetic() == wide &&
            profile.metadata().arithmetic_id == "F02/longdouble-cpu/v1",
        "moved cache retains arithmetic identity");
  std::printf("{\"suite\":\"independent_longdouble_backend\",\"checks\":%d,"
              "\"rational_cast_error\":%.17g,\"passed\":true}\n",
              checks, x.output_rounding_error_inf);
}
