#include "../src/payload_accounting.hpp"
#include "irred/numerics.hpp"
#include "irred/statistics.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred;
int main() {
  unsigned checks = 0;
  auto check = [&](bool b, const char *s) {
    ++checks;
    if (!b)
      throw std::runtime_error(s);
  };
  try {
    detail::PayloadAccounting a(17);
    a.embedded(41, 13);
    check(a.result() == 45, "embedded header once");
    detail::PayloadAccounting b(17);
    b.embedded(12, 13);
    check(!b.result(), "invalid embedded child");
    detail::PayloadAccounting d(17);
    d.embedded({}, 13);
    check(!d.result(), "child overflow propagates");
    detail::PayloadAccounting e(SIZE_MAX - 2);
    e.add(3, 1);
    check(!e.result(), "addition overflow");
    detail::PayloadAccounting f(0);
    f.add(SIZE_MAX, 2);
    check(!f.result(), "multiplication overflow");
    std::vector<double> reserved;
    reserved.reserve(73);
    detail::PayloadAccounting g(11);
    g.vector(reserved);
    check(g.result() == 11 + reserved.capacity() * sizeof(double),
          "empty reserved vector storage");
    std::string text = "x";
    text.reserve(1007);
    detail::PayloadAccounting h(0);
    h.string(text);
    check(h.result() == text.capacity() + 1, "string spare capacity");
    double C[]{4, 1, 1, 9};
    for (auto arithmetic : {numerics::Arithmetic::binary64_legacy_v1,
                            numerics::Arithmetic::longdouble_cpu_v1}) {
      auto factor = numerics::cholesky(C, 2, 4, arithmetic);
      check(factor.status() == numerics::Status::ok, "factor");
      size_t storage = arithmetic == numerics::Arithmetic::binary64_legacy_v1
                           ? 8 * sizeof(double)
                           : 4 * (sizeof(double) + sizeof(long double));
      check(factor.retained_payload_bound() == sizeof(factor) + storage,
            "known two retained matrices");
      auto copied = factor;
      check(copied.retained_payload_bound() == factor.retained_payload_bound(),
            "copy owns same payload");
      auto moved = std::move(copied);
      check(moved.retained_payload_bound() == factor.retained_payload_bound(),
            "move retained payload");
      double rhs[]{5, 10};
      auto sol = numerics::solve(moved, rhs, 1e-8);
      check(sol.status == numerics::Status::ok && sol.value[0] == 1 &&
                sol.value[1] == 1,
            "known solution unaffected");
      statistics::Metadata m;
      m.ordered_ids = {"a", std::string(87, 'b')};
      m.measure = "explicit synthetic";
      m.ordering_provenance = "caller axis declaration";
      auto gaussian = statistics::prepare_gaussian(
          C, statistics::MatrixKind::covariance, m, 4, 1e-8, arithmetic);
      auto before = gaussian.retained_payload_bound();
      check(before && *before > sizeof(gaussian), "Gaussian retained state");
      double residual[]{1, 2};
      auto density = gaussian.evaluate(residual, m.ordered_ids, 1e-8);
      check(density.density.status == statistics::DensityStatus::finite &&
                std::abs(density.quadratic - .6) < 1e-15,
            "fixed quadratic unaffected");
      check(before == gaussian.retained_payload_bound(),
            "evaluation no retained growth");
      double ones[]{1, 1};
      auto op =
          std::move(gaussian).prepare_offset_profile(ones, m.ordered_ids, 1e-8);
      auto cache = op.retained_payload_bound();
      check(cache && *cache > sizeof(op), "profile retained state");
      auto prof = op.evaluate(residual, m.ordered_ids, 1e-8);
      check(prof.status == statistics::DensityStatus::finite &&
                std::abs(prof.quadratic - 1. / 11) < 1e-15,
            "independent profiled contrast");
      check(cache == op.retained_payload_bound(),
            "profile evaluation no retained growth");
    }
    std::cout << "Independent retained payload PASS " << checks << '\n';
  } catch (const std::exception &x) {
    std::cerr << "FAIL " << checks << ": " << x.what() << '\n';
    return 1;
  }
}
