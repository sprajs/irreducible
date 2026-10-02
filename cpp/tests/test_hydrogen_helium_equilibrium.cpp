#include "irred/hydrogen_equilibrium.hpp"
#include "irred/hydrogen_helium_equilibrium.hpp"
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
using W = long double;
using S = irred::numerics::Status;
unsigned checks = 0;
void check(bool ok, const char *what) {
  ++checks;
  if (!ok)
    throw std::runtime_error(what);
}
void close(W value, W expected, W tolerance, const char *what) {
  check(value > 0 && expected > 0 && std::isfinite(value) &&
            std::abs(value / expected - 1) <= tolerance,
        what);
}
std::array<const a::HydrogenHeliumValue *, 6>
groups(const a::HydrogenHeliumRow &r) {
  return {&r.hydrogen_neutral,      &r.hydrogen_ionized,
          &r.helium_neutral,        &r.helium_singly_ionized,
          &r.helium_doubly_ionized, &r.electron_density};
}
W q_e(W T) {
  return std::pow(2 * std::numbers::pi_v<W> * a::hydrogen_electron_mass_kg *
                      irred::boltzmann_constant_joule_per_kelvin * T /
                      (irred::planck_constant_joule_second *
                       irred::planck_constant_joule_second),
                  1.5L);
}
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
} // namespace
int main() {
  try {
    static_assert(a::helium_first_ionization_energy_ev == 24.587389011L);
    static_assert(a::helium_second_ionization_energy_ev == 54.4177655282L);
    std::array states{a::HydrogenHeliumState{3000, 1e8, 1e7},
                      a::HydrogenHeliumState{10000, 6.769876143152199e20, 1e19},
                      a::HydrogenHeliumState{30000, 1e20, 8e18},
                      a::HydrogenHeliumState{100000, 1e20, 1e19}};
    auto out = a::evaluate_hydrogen_helium_equilibrium(states);
    check(out.status == S::ok && out.rows.size() == states.size() &&
              out.solves == states.size(),
          "ordered coarse batch");
    std::size_t iterations = 0, evaluations = 0;
    for (std::size_t i = 0; i < states.size(); ++i) {
      const auto &r = out.rows[i];
      const auto &s = states[i];
      check(r.admission_status == S::ok && r.solves == 1 &&
                r.root_iterations > 0 && r.root_iterations <= 256 &&
                r.charge_evaluations == r.root_iterations + 1,
            "actual bounded row work");
      iterations += r.root_iterations;
      evaluations += r.charge_evaluations;
      for (unsigned j = 0; j < 6; ++j) {
        const auto &v = *groups(r)[j];
        check(v.status == S::ok && v.value && *v.value > 0 &&
                  (j == 5 || *v.value <= 1),
              "independent positive requested group");
        check(v.relative_arithmetic_estimate > 0 &&
                  v.relative_arithmetic_estimate <= 1e-13,
              "per-group diagnostic admission");
      }
      const W h0 = *r.hydrogen_neutral.value, h1 = *r.hydrogen_ionized.value,
              f0 = *r.helium_neutral.value, f1 = *r.helium_singly_ionized.value,
              f2 = *r.helium_doubly_ionized.value,
              e = *r.electron_density.value;
      close(h0 + h1, 1, 2e-15L, "hydrogen nuclei conservation");
      close(f0 + f1 + f2, 1, 2e-15L, "helium nuclei conservation");
      close(e,
            s.hydrogen_nuclei_per_cubic_metre * h1 +
                s.helium_nuclei_per_cubic_metre * (f1 + 2 * f2),
            2e-12L, "shared electron neutrality");
      const W kT =
          irred::boltzmann_constant_joule_per_kelvin * s.temperature_kelvin;
      // Chemical-potential identities witness the three ground degeneracy
      // factors, independent of production root iteration and normalization.
      const W H = q_e(s.temperature_kelvin) *
                  std::exp(-a::hydrogen_ionization_energy_ev *
                           irred::electron_volt_joule / kT),
              He1 = q_e(s.temperature_kelvin) *
                    std::exp(-a::helium_first_ionization_energy_ev *
                             irred::electron_volt_joule / kT),
              He2 = q_e(s.temperature_kelvin) *
                    std::exp(-a::helium_second_ionization_energy_ev *
                             irred::electron_volt_joule / kT);
      close(e * h1 / h0, H, 2e-12L, "degeneracy witness H prefactor1");
      close(e * f1 / f0, 4 * He1, 2e-12L, "degeneracy witness HeI prefactor4");
      close(e * f2 / f1, He2, 2e-12L, "degeneracy witness HeII prefactor1");
      check(std::abs(e * f1 / f0 / He1 - 1) > 2.9L,
            "wrong HeI unit-prefactor refused by control");
      check(bits(r.source.temperature_kelvin) == bits(s.temperature_kelvin) &&
                bits(r.source.hydrogen_nuclei_per_cubic_metre) ==
                    bits(s.hydrogen_nuclei_per_cubic_metre) &&
                bits(r.source.helium_nuclei_per_cubic_metre) ==
                    bits(s.helium_nuclei_per_cubic_metre),
            "source bits/order");
    }
    check(out.root_iterations == iterations &&
              out.charge_evaluations == evaluations,
          "cumulative actual work accounting");
    // Fourth factor witness: helium contributes to the SAME electron owner;
    // independently solving pure hydrogen would violate mixture equilibrium.
    std::array coupled{a::HydrogenHeliumState{30000, 1e19, 1e21}};
    auto mixture = a::evaluate_hydrogen_helium_equilibrium(coupled);
    std::array pure{a::HydrogenState{30000, 1e19}};
    auto hydrogen = a::evaluate_hydrogen_equilibrium(pure);
    check(mixture.rows[0].hydrogen_ionized.value &&
              hydrogen.rows[0].ionized.value &&
              *mixture.rows[0].hydrogen_ionized.value <
                  *hydrogen.rows[0].ionized.value &&
              *mixture.rows[0].electron_density.value >
                  10 * 1e19 * *hydrogen.rows[0].ionized.value,
          "degeneracy/shared-neutrality fourth witness");
    const auto owned = out.rows[0].source;
    states[0] = {9000, 123, 456};
    check(bits(out.rows[0].source.temperature_kelvin) ==
                  bits(owned.temperature_kelvin) &&
              bits(out.rows[0].source.helium_nuclei_per_cubic_metre) ==
                  bits(owned.helium_nuclei_per_cubic_metre),
          "source ownership after caller mutation");
    for (unsigned mask = 1; mask < 64; ++mask) {
      a::HydrogenHeliumPolicy p;
      p.requested_outputs = mask;
      auto result = a::evaluate_hydrogen_helium_equilibrium(
          std::span(states).subspan(1, 1), p);
      check(result.rows[0].admission_status == S::ok,
            "mixed mask root available");
      for (unsigned j = 0; j < 6; ++j) {
        const auto &v = *groups(result.rows[0])[j];
        check(mask & (1u << j)
                  ? v.status == S::ok && v.value.has_value()
                  : !v.value && v.relative_arithmetic_estimate == 0,
              "requested-only group omission");
      }
    }
    std::array tails{a::HydrogenHeliumState{1, 1, 1},
                     a::HydrogenHeliumState{100000, 1e-280, 1e-280},
                     a::HydrogenHeliumState{100000, 1e-290, 1e-290}};
    auto tail = a::evaluate_hydrogen_helium_equilibrium(tails);
    check(tail.rows[0].admission_status == S::ok &&
              tail.rows[0].hydrogen_neutral.value == 1 &&
              tail.rows[0].helium_neutral.value == 1 &&
              !tail.rows[0].electron_density.value &&
              tail.rows[0].electron_density.status ==
                  S::conditioning_budget_exceeded &&
              !tail.rows[0].hydrogen_ionized.value &&
              !tail.rows[0].helium_singly_ionized.value &&
              !tail.rows[0].helium_doubly_ionized.value,
          "tiny charged outputs refused while neutrals available");
    check(tail.rows[1].hydrogen_neutral.value &&
              *tail.rows[1].hydrogen_neutral.value > 0 &&
              tail.rows[1].hydrogen_ionized.value == 1 &&
              tail.rows[1].helium_doubly_ionized.value == 1 &&
              !tail.rows[1].helium_neutral.value,
          "separate subnormal and dominant partners");
    check(tail.rows[2].hydrogen_neutral.status ==
                  S::conditioning_budget_exceeded &&
              !tail.rows[2].hydrogen_neutral.value &&
              tail.rows[2].hydrogen_ionized.value == 1 &&
              tail.rows[2].electron_density.value,
          "insufficient positive cast precision refused separately");
    a::HydrogenHeliumPolicy p;
    p.requested_outputs =
        a::hydrogen_neutral_fraction | a::helium_neutral_fraction;
    auto one =
        a::evaluate_hydrogen_helium_equilibrium(std::span(tails).first(1), p);
    check(one.rows[0].hydrogen_neutral.value == 1 &&
              one.rows[0].helium_neutral.value == 1 &&
              !one.rows[0].electron_density.value &&
              one.rows[0].electron_density.relative_arithmetic_estimate == 0,
          "neutral-only groups do not materialize electron density");
    std::array unrepresented_e{a::HydrogenHeliumState{1, 1e7, 1e6}};
    one = a::evaluate_hydrogen_helium_equilibrium(unrepresented_e, p);
    check(one.rows[0].admission_status == S::ok &&
              one.rows[0].hydrogen_neutral.value == 1 &&
              one.rows[0].helium_neutral.value == 1 &&
              !one.rows[0].electron_density.value,
          "log charge derivative works below wide electron representation");
    const double nan = std::bit_cast<double>(std::uint64_t{0x7ff8000000000042});
    std::array<a::HydrogenHeliumState, 9> bad{{{0, 1, 1},
                                               {100001, 1, 1},
                                               {10000, -0., 1},
                                               {10000, 1, -0.},
                                               {nan, 1, 1},
                                               {10000, INFINITY, 1},
                                               {10000, 1, nan},
                                               {10000, 1e30, 1},
                                               {3000, 1e8, 1e7}}};
    auto invalid = a::evaluate_hydrogen_helium_equilibrium(bad);
    check(invalid.status == S::ok && invalid.solves == 1 &&
              invalid.rows.size() == bad.size(),
          "mixed status order");
    for (unsigned i : {0u, 1u, 2u, 3u, 7u})
      check(invalid.rows[i].admission_status == S::outside_domain &&
                invalid.rows[i].solves == 0,
            "domain refused before solve");
    for (unsigned i : {4u, 5u, 6u})
      check(invalid.rows[i].admission_status == S::nonfinite_input &&
                invalid.rows[i].solves == 0,
            "nonfinite cause retained");
    check(bits(invalid.rows[2].source.hydrogen_nuclei_per_cubic_metre) ==
                  bits(-0.) &&
              bits(invalid.rows[6].source.helium_nuclei_per_cubic_metre) ==
                  bits(nan),
          "failed input signed-zero/NaN bits retained");
    W quantum = q_e(10000);
    std::array boundary{a::HydrogenHeliumState{10000, double(.00045L * quantum),
                                               double(.000225L * quantum)},
                        a::HydrogenHeliumState{10000, double(.00055L * quantum),
                                               double(.000275L * quantum)}};
    one = a::evaluate_hydrogen_helium_equilibrium(boundary);
    check(one.rows[0].admission_status == S::ok &&
              one.rows[1].admission_status == S::outside_domain &&
              one.rows[1].solves == 0,
          "conservative combined degeneracy gate sides");
    if (std::numeric_limits<W>::digits == 64) {
      std::array edge{a::HydrogenHeliumState{
          0x1.38e8p+13, 0x1.0029c8c569159p+80, 0x1.0029c8c569159p+79}};
      auto uncertain = a::evaluate_hydrogen_helium_equilibrium(edge);
      check(uncertain.rows[0].admission_status ==
                    S::conditioning_budget_exceeded &&
                uncertain.rows[0].solves == 0,
            "unresolved combined degeneracy threshold refused");
    }
    std::array repeated{states[1], states[1], bad[0]};
    p = {};
    p.maximum_solves = 1;
    one = a::evaluate_hydrogen_helium_equilibrium(repeated, p);
    check(one.solves == 1 && one.rows[0].electron_density.value &&
              one.rows[1].admission_status == S::work_limit &&
              one.rows[2].admission_status == S::outside_domain,
          "cumulative solve cap and invalid cause precedence");
    auto baseline_work =
        a::evaluate_hydrogen_helium_equilibrium(std::span(repeated).first(1));
    p = {};
    p.maximum_charge_evaluations = baseline_work.charge_evaluations + 1;
    one = a::evaluate_hydrogen_helium_equilibrium(repeated, p);
    check(one.rows[0].electron_density.value &&
              one.rows[1].admission_status == S::work_limit &&
              one.charge_evaluations == p.maximum_charge_evaluations &&
              one.rows[1].charge_evaluations == 1 &&
              one.rows[2].admission_status == S::outside_domain,
          "cumulative charged evaluation cap preserves prior output");
    for (std::size_t iterations_cap : {0u, 1u}) {
      p = {};
      p.maximum_root_iterations = iterations_cap;
      one = a::evaluate_hydrogen_helium_equilibrium(
          std::span(repeated).first(1), p);
      check(one.rows[0].admission_status == S::work_limit &&
                !one.rows[0].electron_density.value &&
                one.root_iterations <= iterations_cap,
            "per-row iterations capped without pretending convergence");
    }
    p = {};
    p.maximum_charge_evaluations = 0;
    one = a::evaluate_hydrogen_helium_equilibrium(states, p);
    check(one.charge_evaluations == 0 &&
              one.rows[0].admission_status == S::work_limit,
          "zero evaluation cap");
    p = {};
    p.maximum_rows = 0;
    check(a::evaluate_hydrogen_helium_equilibrium(states, p).status ==
              S::work_limit,
          "configured row cap");
    p = {};
    p.maximum_native_bytes =
        sizeof(a::HydrogenHeliumBatch) + sizeof(a::HydrogenHeliumRow);
    check(a::evaluate_hydrogen_helium_equilibrium(std::span(states).first(1), p)
                  .status == S::ok,
          "exact payload admitted");
    --p.maximum_native_bytes;
    check(a::evaluate_hydrogen_helium_equilibrium(std::span(states).first(1), p)
                  .status == S::work_limit,
          "payload short by one refused");
    for (unsigned mask : {0u, 64u}) {
      p = {};
      p.requested_outputs = mask;
      check(a::evaluate_hydrogen_helium_equilibrium(states, p).status ==
                S::invalid_input,
            "invalid mask");
    }
    for (unsigned field = 0; field < 5; ++field) {
      p = {};
      if (field == 0)
        p.maximum_rows = 65537;
      if (field == 1)
        p.maximum_solves = 65537;
      if (field == 2)
        p.maximum_root_iterations = 257;
      if (field == 3)
        p.maximum_charge_evaluations = 4000001;
      if (field == 4)
        p.maximum_native_bytes = (1ull << 30) + 1;
      check(a::evaluate_hydrogen_helium_equilibrium(states, p).status ==
                S::invalid_input,
            "hard policy cap");
    }
    std::vector<a::HydrogenHeliumState> oversized(65537, {nan, nan, nan});
    fail_at = 0;
    allocations = 0;
    auto rejected = a::evaluate_hydrogen_helium_equilibrium(oversized);
    fail_at = -1;
    check(rejected.status == S::work_limit && rejected.rows.empty() &&
              allocations == 0,
          "bound before scan/allocation");
    check(a::evaluate_hydrogen_helium_equilibrium({}).status == S::ok,
          "empty batch");
    const int rounding = std::fegetround();
    std::fesetround(FE_UPWARD);
    auto unsupported = a::evaluate_hydrogen_helium_equilibrium(states);
    std::fesetround(rounding);
    check(unsupported.status == S::invalid_input && unsupported.rows.empty(),
          "unsupported rounding before equations/allocation");
    const long live_baseline = live;
    fail_at = 1000;
    allocations = 0;
    {
      auto result = a::evaluate_hydrogen_helium_equilibrium(states);
      check(result.status == S::ok, "allocation baseline");
    }
    const long sites = allocations;
    fail_at = -1;
    check(live == live_baseline && sites == 1,
          "one row allocation with RAII cleanup");
    for (long i = 0; i < sites; ++i) {
      fail_at = i;
      allocations = 0;
      bool threw = false;
      try {
        auto result = a::evaluate_hydrogen_helium_equilibrium(states);
      } catch (const std::bad_alloc &) {
        threw = true;
      }
      fail_at = -1;
      check(threw && live == live_baseline,
            "allocation exception propagated without leak");
    }
    std::cout << "PASS " << checks
              << " H-He equilibrium owner controls allocation_sites=" << sites
              << " named_charge_evaluations=" << evaluations << '\n';
  } catch (const std::exception &e) {
    fail_at = -1;
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
