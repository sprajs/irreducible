// Original independent logit-equation bisection, plus frozen 110-decimal-digit
// Decimal logit controls. No production quadratic/quantum-density helpers.
// Atomic/SI assets are shared ancestry; CPU libm bisection alone is not an
// independent arithmetic reference. No external implementation code copied.
#include "irred/hydrogen_equilibrium.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
namespace {
namespace a = irred::atomic;
using W = long double;
using S = irred::numerics::Status;
constexpr W me = 9.1093837139e-31L, kb = 1.380649e-23L, h = 6.62607015e-34L,
            ev = 1.602176634e-19L, chi = 13.598434599702L * ev;
unsigned checks = 0;
void check(bool ok, const char *what) {
  ++checks;
  if (!ok)
    throw std::runtime_error(what);
}
void near(const a::HydrogenFraction &v, W expected, const char *what) {
  check(v.status == S::ok && v.value && *v.value > 0, what);
  if (std::abs(W(*v.value) - expected) > 2e-12L * expected)
    std::cerr << what << " actual=" << *v.value
              << " expected=" << double(expected) << '\n';
  check(std::abs(W(*v.value) - expected) <= 2e-12L * expected, what);
}
std::array<W, 2> oracle(const a::HydrogenState &s, unsigned iterations) {
  const W t = s.temperature_kelvin, n = s.hydrogen_nuclei_per_cubic_metre,
          prefactor =
              std::pow(2 * std::numbers::pi_v<W> * me * kb * t / (h * h), 1.5L),
          target = std::log(prefactor) - chi / (kb * t) - std::log(n);
  // u=ln(x/y). x=logistic(u), y=logistic(-u); solve the original
  // 2 ln x-ln y=ln K monotonically, independent of the quadratic branches.
  W lo = -2000, hi = 2000;
  for (unsigned i = 0; i < iterations; ++i) {
    const W mid = lo + (hi - lo) / 2;
    if (mid == lo || mid == hi)
      break;
    const W x = 1 / (1 + std::exp(-mid)), y = 1 / (1 + std::exp(mid));
    if (2 * std::log(x) - std::log(y) < target)
      lo = mid;
    else
      hi = mid;
  }
  const W u = lo + (hi - lo) / 2;
  return {1 / (1 + std::exp(-u)), 1 / (1 + std::exp(u))};
}
struct Fact {
  a::HydrogenState source;
  W x, y;
};
// Inputs are exact binary64 values. References: Decimal precision110, original
// 800-step logit bisection, independently pinned 110-digit pi and quantum
// density power; repeated precision150 agrees within 1e-100 relative for these
// separately represented fractions. Fraction rounded1 is never a zero partner.
constexpr std::array facts{
    Fact{{1000, 1e6},
         4.73023852832976752438021452630805498865e-25L,
         0.99999999999999999999999952697614716702325L},
    Fact{{3000, 1e8},
         0.007506758504272591528363283130911396733257L,
         0.992493241495727408471636716869088603266743L},
    Fact{{10000, 1e20},
         0.807408862965693455667027652246600477303731L,
         0.192591137034306544332972347753399522696269L},
    Fact{{100000, 1e-280},
         1.L,
         6.34557246361083667034822199409292334484001e-309L},
    Fact{{100000, 1e-270},
         1.L,
         6.34557246361083720633096883401955245030167e-299L}};
} // namespace
int main() {
  try {
    std::array<a::HydrogenState, facts.size()> states{};
    for (unsigned i = 0; i < facts.size(); ++i)
      states[i] = facts[i].source;
    auto result = a::evaluate_hydrogen_equilibrium(states);
    check(result.status == S::ok && result.solves == states.size(),
          "peer batch");
    for (unsigned i = 0; i < facts.size(); ++i) {
      check(result.rows[i].admission_status == S::ok,
            "peer physical admission");
      near(result.rows[i].ionized, facts[i].x,
           "110-digit independent ionized fact");
      near(result.rows[i].neutral, facts[i].y,
           "110-digit independent neutral fact");
      const auto coarse = oracle(states[i], 128), fine = oracle(states[i], 256);
      for (unsigned j = 0; j < 2; ++j)
        check(fine[j] > 0 && std::abs(coarse[j] - fine[j]) <= 2e-13L * fine[j],
              "independent logit refinement allocation");
      near(result.rows[i].ionized, fine[0], "independent log equation ionized");
      near(result.rows[i].neutral, fine[1], "independent log equation neutral");
    }
    // Monotonic temperature/density dependence at admitted fixed controls;
    // ordering/permutation is not extra independent arithmetic evidence.
    std::array grid{a::HydrogenState{3000, 1e8}, a::HydrogenState{4000, 1e8},
                    a::HydrogenState{5000, 1e8}, a::HydrogenState{4000, 1e9}};
    auto out = a::evaluate_hydrogen_equilibrium(grid);
    check(*out.rows[0].ionized.value < *out.rows[1].ionized.value &&
              *out.rows[1].ionized.value < *out.rows[2].ionized.value &&
              *out.rows[3].ionized.value < *out.rows[1].ionized.value,
          "equilibrium monotonic temperature/density");
    std::array reordered{grid[3], grid[1], grid[0], grid[2]};
    auto permuted = a::evaluate_hydrogen_equilibrium(reordered);
    for (unsigned i = 0; i < 4; ++i) {
      const unsigned original = std::array{3u, 1u, 0u, 2u}[i];
      check(permuted.rows[i].ionized.value ==
                    out.rows[original].ionized.value &&
                permuted.rows[i].neutral.value ==
                    out.rows[original].neutral.value,
            "source row permutation equivariance");
    }
    // Asymptotic K<<1: x/sqrt(K)->1; K>>1: K*y->1.
    const W t = 1000, n = 1e6,
            k = std::pow(2 * std::numbers::pi_v<W> * me * kb * t / (h * h),
                         1.5L) *
                std::exp(-chi / (kb * t)) / n;
    check(std::abs(W(*result.rows[0].ionized.value) / std::sqrt(k) - 1) <
              2e-12L,
          "neutral-gas leading limit");
    const auto high = facts[4].source;
    const W highk = std::pow(2 * std::numbers::pi_v<W> * me * kb *
                                 high.temperature_kelvin / (h * h),
                             1.5L) *
                    std::exp(-chi / (kb * high.temperature_kelvin)) /
                    high.hydrogen_nuclei_per_cubic_metre;
    check(std::abs(W(*result.rows[4].neutral.value) * highk - 1) < 2e-12L,
          "ionized-gas leading neutral limit");
    std::cout << "PASS " << checks
              << " hydrogen equilibrium independent peer controls logit128/256 "
                 "Decimal110/150\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
