#include "irred/hydrogen_equilibrium.hpp"
#include "irred/quantities.hpp"
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <numbers>
#include <stdexcept>
static long fail_at = -1, allocations = 0, live = 0;
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(std::size_t n) {
  if (fail_at >= 0 && allocations++ == fail_at)
    throw std::bad_alloc();
  if (void *p = std::malloc(n ? n : 1)) {
    ++live;
    return p;
  }
  throw std::bad_alloc();
}
NOINLINE void *operator new[](std::size_t n) { return ::operator new(n); }
NOINLINE void operator delete(void *p) noexcept {
  if (p) {
    --live;
    std::free(p);
  }
}
NOINLINE void operator delete[](void *p) noexcept { ::operator delete(p); }
NOINLINE void operator delete(void *p, std::size_t) noexcept {
  ::operator delete(p);
}
NOINLINE void operator delete[](void *p, std::size_t) noexcept {
  ::operator delete(p);
}
namespace {
namespace a = irred::atomic;
using S = irred::numerics::Status;
unsigned checks = 0;
void check(bool ok, const char *what) {
  ++checks;
  if (!ok)
    throw std::runtime_error(what);
}
void near(const a::HydrogenFraction &v, long double expected,
          const char *what) {
  check(v.status == S::ok && v.value && *v.value > 0 && *v.value <= 1, what);
  check(std::abs(static_cast<long double>(*v.value) - expected) <=
            2e-12L * std::abs(expected),
        what);
  check(v.relative_arithmetic_estimate > 0 &&
            v.relative_arithmetic_estimate <= 1e-13,
        "fraction arithmetic diagnostic");
}
} // namespace
int main() {
  try {
    static_assert(irred::boltzmann_constant_joule_per_kelvin == 1.380649e-23L);
    static_assert(irred::electron_volt_joule == 1.602176634e-19L);
    // Exact rational x controls: density derived at 110 decimal digits, then
    // frozen as binary64. Rounding shifts x by <1e-16 relative. Public peer
    // retains separate exact-input high-precision regimes, independent of this
    // rational algebra.
    const std::array densities{2.0309628429456599e21, 6.769876143152199e20,
                               2.5387035536820748e20, 4.178935890834691e19};
    const std::array<long double, 4> exact{1.L / 3, .5L, 2.L / 3, .9L};
    std::array<a::HydrogenState, 4> states{};
    for (unsigned i = 0; i < states.size(); ++i)
      states[i] = {10000, densities[i]};
    auto out = a::evaluate_hydrogen_equilibrium(states);
    check(out.status == S::ok && out.solves == 4 && out.rows.size() == 4,
          "coarse rational batch");
    for (unsigned i = 0; i < states.size(); ++i) {
      check(out.rows[i].admission_status == S::ok && out.rows[i].solves == 1,
            "rational row admission");
      near(out.rows[i].ionized, exact[i], "exact rational ionized fraction");
      near(out.rows[i].neutral, 1 - exact[i],
           "exact rational neutral fraction");
      check(std::bit_cast<std::uint64_t>(
                out.rows[i].source.hydrogen_nuclei_per_cubic_metre) ==
                std::bit_cast<std::uint64_t>(
                    states[i].hydrogen_nuclei_per_cubic_metre),
            "exact density bits/order");
    }
    const auto retained = out.rows[0].source;
    states[0] = {20000, 1};
    check(out.rows[0].source.temperature_kelvin ==
                  retained.temperature_kelvin &&
              out.rows[0].source.hydrogen_nuclei_per_cubic_metre ==
                  retained.hydrogen_nuclei_per_cubic_metre,
          "result owns source values after caller mutation");
    a::HydrogenPolicy p;
    p.requested_outputs = a::ionized_fraction;
    auto one =
        a::evaluate_hydrogen_equilibrium(std::span(states).subspan(1, 1), p);
    check(one.rows[0].ionized.value && !one.rows[0].neutral.value &&
              one.rows[0].neutral.relative_arithmetic_estimate == 0,
          "output omission");
    // Independent representations avoid both 1-x cancellation and fake zero.
    std::array tails{
        a::HydrogenState{1000, 1e6}, a::HydrogenState{100000, 1e-280},
        a::HydrogenState{1, 1}, a::HydrogenState{100000, 1e-290},
        a::HydrogenState{100000, std::numeric_limits<double>::denorm_min()}};
    auto thin = a::evaluate_hydrogen_equilibrium(tails);
    check(thin.rows[0].ionized.value && *thin.rows[0].ionized.value > 0 &&
              thin.rows[0].neutral.value == 1,
          "tiny ionized and rounded-one neutral coexist");
    check(thin.rows[1].neutral.value && *thin.rows[1].neutral.value > 0 &&
              thin.rows[1].ionized.value == 1,
          "subnormal neutral retained separately");
    check(thin.rows[2].admission_status == S::ok &&
              !thin.rows[2].ionized.value &&
              thin.rows[2].ionized.status == S::conditioning_budget_exceeded &&
              thin.rows[2].neutral.value == 1,
          "positive ionization underflow never zero");
    for (unsigned i : {3u, 4u})
      check(!thin.rows[i].neutral.value &&
                thin.rows[i].neutral.status ==
                    S::conditioning_budget_exceeded &&
                thin.rows[i].ionized.value == 1,
            "tiny neutral cast loss/underflow refusal");
    p.requested_outputs = a::neutral_fraction;
    one = a::evaluate_hydrogen_equilibrium(std::span(tails).subspan(2, 1), p);
    check(one.rows[0].neutral.value == 1 && !one.rows[0].ionized.value,
          "rounded neutral available when tiny partner unrequested");
    // Sufficient total-density guard comes before fraction execution. Source
    // exactness and each invalid row cause survive a mixed batch.
    const double nan = std::bit_cast<double>(std::uint64_t{0x7ff8000000000042});
    std::array bad{a::HydrogenState{0, 1},
                   a::HydrogenState{1e5 + 1, 1},
                   a::HydrogenState{10000, -0.},
                   a::HydrogenState{nan, 1},
                   a::HydrogenState{10000, INFINITY},
                   a::HydrogenState{10000, 1e30},
                   a::HydrogenState{3000, 1e8}};
    out = a::evaluate_hydrogen_equilibrium(bad);
    check(out.status == S::ok && out.solves == 1 &&
              out.rows.size() == bad.size(),
          "mixed batch preserves good/invalid row order");
    for (unsigned i : {0u, 1u, 2u, 5u})
      check(out.rows[i].admission_status == S::outside_domain &&
                out.rows[i].solves == 0,
            "physical domain refused before solve");
    for (unsigned i : {3u, 4u})
      check(out.rows[i].admission_status == S::nonfinite_input,
            "nonfinite row tagged");
    check(std::bit_cast<std::uint64_t>(
              out.rows[2].source.hydrogen_nuclei_per_cubic_metre) ==
                  std::bit_cast<std::uint64_t>(-0.) &&
              std::bit_cast<std::uint64_t>(
                  out.rows[3].source.temperature_kelvin) ==
                  std::bit_cast<std::uint64_t>(nan),
          "failed source signed zero and NaN payload bits");
    const long double base = 2 * std::numbers::pi_v<long double> *
                             a::hydrogen_electron_mass_kg *
                             irred::boltzmann_constant_joule_per_kelvin *
                             10000 /
                             (irred::planck_constant_joule_second *
                              irred::planck_constant_joule_second),
                      nq = base * std::sqrt(base);
    std::array boundary{a::HydrogenState{10000, double(.0009L * nq)},
                        a::HydrogenState{10000, double(.0011L * nq)}};
    one = a::evaluate_hydrogen_equilibrium(boundary);
    check(one.rows[0].admission_status == S::ok &&
              one.rows[1].admission_status == S::outside_domain &&
              one.rows[1].solves == 0,
          "nondegeneracy guard sides");
    if (std::numeric_limits<long double>::digits == 64) {
      // Frozen binary64 input within the 64-bit-wide density admission margin;
      // cannot label a scientifically unresolved boundary as admitted/clipped.
      std::array unresolved{
          a::HydrogenState{0x1.38e8p+13, 0x1.0029c8c569159p+81}};
      const auto edge = a::evaluate_hydrogen_equilibrium(unresolved);
      check(edge.rows[0].admission_status == S::conditioning_budget_exceeded &&
                edge.rows[0].solves == 0 && !edge.rows[0].ionized.value,
            "unresolved nondegeneracy boundary refusal");
    }
    p = {};
    p.maximum_solves = 1;
    one = a::evaluate_hydrogen_equilibrium(boundary, p);
    check(one.solves == 1 && one.rows[1].admission_status == S::outside_domain,
          "physical invalid cause precedes exhausted work");
    states[0] = states[1];
    out = a::evaluate_hydrogen_equilibrium(states, p);
    check(out.solves == 1 && out.rows[0].ionized.value &&
              out.rows[1].admission_status == S::work_limit &&
              !out.rows[1].ionized.value,
          "global solves preserve prior row");
    p = {};
    p.maximum_rows = 0;
    check(a::evaluate_hydrogen_equilibrium(states, p).status == S::work_limit,
          "row quota");
    p = {};
    p.maximum_native_bytes = sizeof(a::HydrogenBatch) + sizeof(a::HydrogenRow);
    check(a::evaluate_hydrogen_equilibrium(std::span(states).first(1), p)
                  .status == S::ok,
          "exact result payload threshold");
    --p.maximum_native_bytes;
    check(a::evaluate_hydrogen_equilibrium(std::span(states).first(1), p)
                  .status == S::work_limit,
          "result payload short by one");
    for (unsigned mask : {0u, 4u}) {
      p = {};
      p.requested_outputs = mask;
      check(a::evaluate_hydrogen_equilibrium(states, p).status ==
                S::invalid_input,
            "invalid output mask");
    }
    p = {};
    p.maximum_rows = 65537;
    check(a::evaluate_hydrogen_equilibrium(states, p).status ==
              S::invalid_input,
          "hard policy row maximum");
    p = {};
    p.maximum_solves = 65537;
    check(a::evaluate_hydrogen_equilibrium(states, p).status ==
              S::invalid_input,
          "hard solve policy maximum");
    p = {};
    p.maximum_native_bytes = (1ull << 30) + 1;
    check(a::evaluate_hydrogen_equilibrium(states, p).status ==
              S::invalid_input,
          "hard byte policy maximum");
    std::vector<a::HydrogenState> oversize(65537, {nan, nan});
    allocations = 0;
    fail_at = 0;
    auto rejected = a::evaluate_hydrogen_equilibrium(oversize);
    fail_at = -1;
    check(rejected.status == S::work_limit && allocations == 0,
          "bounds before scientific scan/allocation");
    check(a::evaluate_hydrogen_equilibrium({}).status == S::ok, "empty batch");
    const int rounding = std::fegetround();
    std::fesetround(FE_UPWARD);
    auto unsupported = a::evaluate_hydrogen_equilibrium(states);
    std::fesetround(rounding);
    check(unsupported.status == S::invalid_input && unsupported.rows.empty(),
          "unsupported rounding before allocation/equations");
    const long baseline = live;
    allocations = 0;
    fail_at = 1000;
    {
      auto result = a::evaluate_hydrogen_equilibrium(states);
      check(result.status == S::ok, "allocation baseline");
    }
    const long sites = allocations;
    fail_at = -1;
    check(live == baseline && sites == 1, "one coarse row allocation cleanup");
    for (long i = 0; i < sites; ++i) {
      allocations = 0;
      fail_at = i;
      bool threw = false;
      try {
        auto result = a::evaluate_hydrogen_equilibrium(states);
      } catch (const std::bad_alloc &) {
        threw = true;
      }
      fail_at = -1;
      check(threw && live == baseline, "allocation failure RAII cleanup");
    }
    std::cout << "PASS " << checks
              << " hydrogen equilibrium owner controls allocation_sites="
              << sites << '\n';
    return 0;
  } catch (const std::exception &e) {
    fail_at = -1;
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
