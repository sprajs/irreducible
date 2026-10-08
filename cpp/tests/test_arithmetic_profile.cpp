// Analytic dyadic SPD solve and explicit unsupported-wide refusal. This test
// checks the actual numerical owner; machine properties alone are not a pass.
#include "irred/numerics.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string_view>

namespace {
void need(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL " << message << '\n';
    std::exit(1);
  }
}
}

int main(int argc, char **argv) {
  using namespace irred::numerics;
  using D = std::numeric_limits<double>;
  using W = std::numeric_limits<long double>;
  if (argc == 2 && std::string_view(argv[1]) == "--expect-binary64-long-double")
    need(W::digits == D::digits && W::max_exponent == D::max_exponent &&
             W::min_exponent == D::min_exponent,
         "declared Apple Silicon long-double representation");
  else
    need(argc == 1, "unknown arithmetic-profile argument");
  need(D::radix == 2 && D::digits == 53 && D::is_iec559,
       "binary64 representation");
  need(std::fegetround() == FE_TONEAREST, "nearest rounding at entry");
  const std::array<double, 4> covariance{4, 2, 2, 5};
  const std::array<double, 2> rhs{0, -4};
  const auto binary = cholesky(covariance, 2, 4);
  need(binary.status() == Status::ok &&
           binary.arithmetic() == Arithmetic::binary64_legacy_v1,
       "binary64 SPD preparation");
  const auto solve_check = [&](const Factorization &factor) {
    const auto result = solve(factor, rhs, 1e-10);
    need(result.status == Status::ok && result.value.size() == 2,
         "analytic SPD solve admitted");
    need(std::abs(result.value[0] - .5) <= 1e-14 &&
             std::abs(result.value[1] + 1) <= 1e-14,
         "independent exact solution (1/2,-1)");
  };
  solve_check(binary);
  const bool wide_supported = W::radix == 2 && W::digits >= 64 &&
      W::min_exponent < D::min_exponent && W::max_exponent > D::max_exponent;
  const auto wide = cholesky(covariance, 2, 4, Arithmetic::longdouble_cpu_v1);
  need(wide.arithmetic() == Arithmetic::longdouble_cpu_v1,
       "requested arithmetic identity retained");
  if (wide_supported) {
    need(wide.status() == Status::ok, "supported wide SPD preparation");
    solve_check(wide);
    need(std::fesetround(FE_UPWARD) == 0, "directed rounding available");
    const auto directed = cholesky(covariance, 2, 4, Arithmetic::longdouble_cpu_v1);
    need(std::fesetround(FE_TONEAREST) == 0, "restore nearest rounding");
    need(directed.status() == Status::outside_domain,
         "unsupported rounding refuses wide arithmetic");
  } else {
    need(wide.status() == Status::outside_domain && wide.size() == 0,
         "unsupported wide representation refuses without factor payload");
    need(solve(wide, rhs, 1e-10).status != Status::ok,
         "refused wide owner cannot publish a solve");
  }
  std::cout << "{\"suite\":\"arithmetic_profile\",\"double_digits\":"
            << D::digits << ",\"long_double_digits\":" << W::digits
            << ",\"long_double_min_exponent\":" << W::min_exponent
            << ",\"long_double_max_exponent\":" << W::max_exponent
            << ",\"long_double_bytes\":" << sizeof(long double)
            << ",\"wide_supported\":" << (wide_supported ? "true" : "false")
            << ",\"passed\":true}\n";
}
