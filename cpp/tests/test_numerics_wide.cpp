// Wide native policy owner tests. Independent exact rational A det18 and
// adjugate [5,-2,1;-2,8,-4;1,-4,11]/18; constructed dyadic r=Cv.
// Predeclared component absolute1e-12; sensitivity policy1e-10. Empirical
// diagnostic, not certified error bound. Test-only LDLT is separate algorithm.
#include "fixtures/ldlt_reference.hpp"
#include "irred/numerics.hpp"
#include "irred/statistics.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
int main() {
  try {
    int count = 0;
    auto check = [&](bool v) {
      ++count;
      if (!v)
        throw std::runtime_error("wide owner contract");
    };
    std::array<double, 9> a{4, 1, 0, 1, 3, 1, 0, 1, 2};
    std::array<double, 3> r{1, -2, 3}, cv{1, -.5, 3}, zero{0, 0, 0};
    auto f =
        numerics::cholesky(a, 3, 9, numerics::Arithmetic::longdouble_cpu_v1);
    check(f.status() == numerics::Status::ok);
    check(f.arithmetic() == numerics::Arithmetic::longdouble_cpu_v1);
    auto x = numerics::solve(f, r, 1e-10);
    check(x.status == numerics::Status::ok);
    auto ref = test_reference::solve(test_reference::factor(a, 3), r);
    for (unsigned i = 0; i < 3; ++i)
      check(std::abs(x.value[i] - ref[i]) < 1e-12L);
    check(x.estimated_forward_sensitivity >= x.output_rounding_error_relative);
    check(x.pre_cast_sensitivity_estimate > 0);
    check(x.post_cast_backward_residual == x.backward_residual);
    check(x.output_rounding_error_inf >= 0);
    auto z = numerics::solve(f, zero, 1e-10);
    check(z.status == numerics::Status::ok &&
          z.output_rounding_error_inf == 0 &&
          z.output_rounding_error_relative == 0);
    auto v = numerics::solve(f, cv, 1e-10);
    check(v.status == numerics::Status::ok);
    for (unsigned i = 0; i < 3; ++i)
      check(std::abs(v.value[i] - std::array<double, 3>{.5, -1, 2}[i]) < 1e-12);
    check(numerics::solve(f, r, 1e-30).status ==
          numerics::Status::conditioning_budget_exceeded);
    check(numerics::cholesky(a, 3, 9, static_cast<numerics::Arithmetic>(999))
              .status() == numerics::Status::invalid_input);
    auto unsupported =
        numerics::cholesky(a, 3, 9, static_cast<numerics::Arithmetic>(999));
    check(std::string(unsupported.arithmetic_id()) ==
          "F02/unsupported-arithmetic");
    std::feclearexcept(FE_ALL_EXCEPT);
    std::feraiseexcept(FE_UNDERFLOW | FE_INEXACT);
    auto flags = std::fetestexcept(FE_ALL_EXCEPT);
    auto preserved = numerics::solve(f, r, 1e-10);
    check(preserved.status == numerics::Status::ok);
    check(std::fetestexcept(FE_ALL_EXCEPT) == flags);
    std::feclearexcept(FE_ALL_EXCEPT);
    auto old = std::fegetround();
    std::fesetround(FE_UPWARD);
    check(numerics::solve(f, r, 1e-10).status ==
          numerics::Status::outside_domain);
    check(numerics::cholesky(a, 3, 9, numerics::Arithmetic::longdouble_cpu_v1)
              .status() == numerics::Status::outside_domain);
    std::fesetround(old);
    std::array<double, 1> one{1},
        tiny{std::numeric_limits<double>::denorm_min()};
    auto scalar =
        numerics::cholesky(one, 1, 1, numerics::Arithmetic::longdouble_cpu_v1);
    check(numerics::solve(scalar, tiny, 1e-10).status ==
          numerics::Status::outside_domain);
    statistics::Metadata m;
    m.ordered_ids = {"a", "b", "c"};
    m.measure = "product dx";
    m.ordering_provenance = "generated explicit";
    auto invalid_gaussian = statistics::prepare_gaussian(
        a, statistics::MatrixKind::covariance, m, 9, 1e-10,
        static_cast<numerics::Arithmetic>(999));
    check(invalid_gaussian.status() ==
              statistics::DensityStatus::invalid_input &&
          invalid_gaussian.metadata().arithmetic_id ==
              "F02/unsupported-arithmetic");
    auto g = statistics::prepare_gaussian(
        a, statistics::MatrixKind::covariance, m, 9, 1e-10,
        numerics::Arithmetic::longdouble_cpu_v1);
    check(g.status() == statistics::DensityStatus::finite);
    std::array<std::size_t, 2> keep{0, 2};
    auto selected = g.marginal(keep, 4, 1e-10);
    check(selected.arithmetic() == numerics::Arithmetic::longdouble_cpu_v1);
    check(selected.metadata().arithmetic_id == "F02/longdouble-cpu/v1");
    std::array<double, 3> response{1, 1, 1};
    auto prior = g.proper_offset(response, m.ordered_ids, 0, 1, "latent", true,
                                 9, 1e-10);
    check(prior.arithmetic() == numerics::Arithmetic::longdouble_cpu_v1);
    auto profile =
        std::move(g).prepare_offset_profile(response, m.ordered_ids, 1e-10);
    check(profile.arithmetic() == numerics::Arithmetic::longdouble_cpu_v1);
    auto score = profile.evaluate(r, m.ordered_ids, 1e-10);
    check(score.status == statistics::DensityStatus::finite &&
          std::abs(score.quadratic - 61.L / 7) < 1e-12L);
    std::printf(
        "{\"suite\":\"native_wide_owner\",\"checks\":%d,\"passed\":true}\n",
        count);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
