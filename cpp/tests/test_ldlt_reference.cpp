// Predeclared reference validation1e-12 abs for exact integer3x3 fixtures.
// A=[4,1,0;1,3,1;0,1,2], det18, adjugate=[5,-2,1;-2,8,-4;1,-4,11].
// r=[1,-2,3] gives x=[2/3,-5/3,7/3], q11; ones Gram7/9, b4/3,
// profile a12/7,q61/7. Constructed v=[1/2,-1,2], Av=[1,-1/2,3]
// uses exact dyadic inputs and has q7. No production factor/solve called.
#include "fixtures/ldlt_reference.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
int main() {
  try {
    using namespace irred::test_reference;
    int checks = 0;
    auto check = [&](bool v) {
      ++checks;
      if (!v)
        throw std::runtime_error("LDLT reference fixture");
    };
    auto near = [&](long double a, long double b) {
      check(std::abs(a - b) <= 1e-12L);
    };
    std::array<double, 9> c{4, 1, 0, 1, 3, 1, 0, 1, 2};
    auto f = factor(c, 3);
    check(f.valid);
    near(f.logdet, std::log(18.L));
    std::array<double, 3> r{1, -2, 3}, one{1, 1, 1};
    auto x = solve(f, r), u = solve(f, one);
    near(x[0], 2.L / 3);
    near(x[1], -5.L / 3);
    near(x[2], 7.L / 3);
    near(quadratic(r, x), 11);
    long double b = 0, g = 0;
    for (unsigned i = 0; i < 3; ++i) {
      b += x[i];
      g += u[i];
    }
    near(b, 4.L / 3);
    near(g, 7.L / 9);
    auto a = b / g;
    near(a, 12.L / 7);
    std::array<double, 3> adjusted{};
    for (unsigned i = 0; i < 3; ++i)
      adjusted[i] = static_cast<double>(r[i] - a);
    near(quadratic(adjusted, solve(f, adjusted)), 61.L / 7);
    std::array<double, 3> cv{1, -.5, 3};
    auto v = solve(f, cv);
    near(v[0], .5);
    near(v[1], -1);
    near(v[2], 2);
    near(quadratic(cv, v), 7);
    auto bad = c;
    bad[1] = 2;
    check(!factor(bad, 3).valid);
    bad = c;
    bad[0] = -1;
    check(!factor(bad, 3).valid);
    check(!factor({}, 0).valid);
    std::printf("{\"suite\":\"independent_LDLT_reference\",\"checks\":%d,"
                "\"passed\":true}\n",
                checks);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
