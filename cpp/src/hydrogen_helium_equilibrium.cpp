#include "irred/hydrogen_helium_equilibrium.hpp"
#include "hydrogen_quantum_density.hpp"
#include "irred/quantities.hpp"
#include "payload_accounting.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::atomic {
namespace {
using W = long double;
using S = numerics::Status;
constexpr std::size_t hard_rows = 65536, hard_evaluations = 4000000,
                      hard_iterations = 256, hard_bytes = 1ull << 30;
constexpr W arithmetic = 64.L * std::numeric_limits<W>::epsilon();
constexpr W admission = 1e-13L;
constexpr W ln2 = 0.693147180559945309417232121458176568L;
std::array<HydrogenHeliumValue *, 6> groups(HydrogenHeliumRow &r) {
  return {&r.hydrogen_neutral,      &r.hydrogen_ionized,
          &r.helium_neutral,        &r.helium_singly_ionized,
          &r.helium_doubly_ionized, &r.electron_density};
}
void fail(HydrogenHeliumRow &r, S status, unsigned mask) {
  r.admission_status = status;
  auto values = groups(r);
  for (unsigned i = 0; i < values.size(); ++i)
    if (mask & (1u << i))
      values[i]->status = status;
}
template <std::size_t N> W log_sum(const std::array<W, N> &terms) {
  W maximum = terms[0];
  for (W v : terms) {
    if (!std::isfinite(v))
      return std::numeric_limits<W>::quiet_NaN();
    maximum = std::max(maximum, v);
  }
  W sum = 0;
  for (W v : terms)
    sum += std::exp(v - maximum);
  return maximum + std::log(sum);
}
template <std::size_t N>
std::array<W, N> normalized_logs(const std::array<W, N> &weights) {
  W maximum = weights[0];
  for (W v : weights)
    maximum = std::max(maximum, v);
  std::array<W, N> shifted{};
  W sum = 0;
  for (std::size_t i = 0; i < N; ++i) {
    shifted[i] = weights[i] - maximum;
    sum += std::exp(shifted[i]);
  }
  // Work with shifted weights even for the dominant stage. Subtracting the
  // unshifted normalization could destroy a small positive log probability.
  const W normalization = std::log(sum);
  for (W &v : shifted)
    v -= normalization;
  return shifted;
}
struct Charge {
  std::array<W, 2> h;
  std::array<W, 3> he;
  W log_density = 0, derivative = 0;
};
struct Equation {
  W log_hydrogen, log_helium;
  std::array<W, 3> log_saha;
  Charge at(W ell) const {
    Charge r{normalized_logs(std::array<W, 2>{0, log_saha[0] - ell}),
             normalized_logs(std::array<W, 3>{
                 0, log_saha[1] - ell, log_saha[1] + log_saha[2] - 2 * ell}),
             0, 0};
    r.log_density =
        log_sum(std::array<W, 3>{log_hydrogen + r.h[1], log_helium + r.he[1],
                                 log_helium + r.he[2] + ln2});
    // dF/dell = 1 + weighted charge variance/C, in positive products.
    r.derivative =
        1 + std::exp(log_hydrogen + r.h[0] + r.h[1] - r.log_density) +
        std::exp(log_helium + r.he[0] + r.he[1] - r.log_density) +
        4 * std::exp(log_helium + r.he[0] + r.he[2] - r.log_density) +
        std::exp(log_helium + r.he[1] + r.he[2] - r.log_density);
    return r;
  }
};
void store(HydrogenHeliumValue &out, W log_value, W relative, bool fraction) {
  const W value = std::exp(log_value);
  const double rounded = static_cast<double>(value);
  if (!std::isfinite(value) || !std::isfinite(rounded) || !(value > 0) ||
      !(rounded > 0) || (fraction && (value > 1 || rounded > 1))) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  relative += arithmetic * (1 + std::abs(log_value)) +
              std::abs(W(rounded) - value) / value;
  if (!std::isfinite(relative) || relative > admission) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out.value = rounded;
  out.relative_arithmetic_estimate = static_cast<double>(relative);
  out.status = S::ok;
}
void solve(HydrogenHeliumRow &row, HydrogenHeliumBatch &out,
           const HydrogenHeliumPolicy &policy, const Equation &eq, W upper,
           const std::array<W, 3> &saha_error) {
  auto evaluate = [&](W ell, Charge &r) {
    if (out.charge_evaluations == policy.maximum_charge_evaluations) {
      fail(row, S::work_limit, policy.requested_outputs);
      return false;
    }
    ++out.charge_evaluations;
    ++row.charge_evaluations;
    r = eq.at(ell);
    if (!std::isfinite(r.log_density) || !std::isfinite(r.derivative) ||
        !(r.derivative >= 1)) {
      fail(row, S::conditioning_budget_exceeded, policy.requested_outputs);
      return false;
    }
    return true;
  };
  if (!policy.maximum_root_iterations) {
    fail(row, S::work_limit, policy.requested_outputs);
    return;
  }
  Charge current;
  if (!evaluate(upper, current))
    return;
  W lower = current.log_density;
  // The lower endpoint is C(e_upper). C decreases with e, so it brackets the
  // unique root without a second owner of the pure-H quadratic equation.
  const W endpoint_margin = arithmetic * (1 + std::abs(upper));
  if (lower > upper + endpoint_margin) {
    fail(row, S::conditioning_budget_exceeded, policy.requested_outputs);
    return;
  }
  lower = std::min(lower, upper);
  W ell = lower + (upper - lower) / 2, residual = 0, root_distance = 0;
  bool converged = false;
  for (std::size_t i = 0; i < policy.maximum_root_iterations; ++i) {
    if (!evaluate(ell, current))
      return;
    ++row.root_iterations;
    ++out.root_iterations;
    residual = ell - current.log_density;
    if (!std::isfinite(residual)) {
      fail(row, S::conditioning_budget_exceeded, policy.requested_outputs);
      return;
    }
    if (residual <= 0)
      lower = ell;
    else
      upper = ell;
    // F' >= 1 bounds root distance by the residual, apart from the explicit
    // empirical evaluation allowance propagated below. Adjacent wide values
    // are also a valid termination with per-output representability refusal.
    root_distance = std::min(std::abs(residual), upper - lower);
    const W target = arithmetic * (1 + std::abs(ell));
    if (std::abs(residual) <= target || upper - lower <= target ||
        std::nextafter(lower, upper) >= upper) {
      converged = true;
      break;
    }
    const W trial = ell - residual / current.derivative;
    const W quarter = (upper - lower) / 4;
    ell = std::isfinite(trial) && trial > lower + quarter &&
                  trial < upper - quarter
              ? trial
              : lower + (upper - lower) / 2;
  }
  if (!converged) {
    fail(row, S::work_limit, policy.requested_outputs);
    return;
  }
  const W h0 = std::exp(current.h[0]), h1 = std::exp(current.h[1]),
          f0 = std::exp(current.he[0]), f1 = std::exp(current.he[1]),
          f2 = std::exp(current.he[2]);
  const W log_mu =
              log_sum(std::array<W, 2>{current.he[1], ln2 + current.he[2]}),
          log_two_minus_mu =
              log_sum(std::array<W, 2>{ln2 + current.he[0], current.he[1]}),
          denominator = current.derivative;
  // Implicit log-electron sensitivities, all positive; use log probabilities
  // so tiny charged stages need not be materialized to qualify neutral groups.
  const W sensitivity_h = std::exp(eq.log_hydrogen + current.h[0] +
                                   current.h[1] - current.log_density) /
                          denominator,
          sensitivity_1 = std::exp(eq.log_helium + current.he[0] + log_mu -
                                   current.log_density) /
                          denominator,
          sensitivity_2 = std::exp(eq.log_helium + current.he[2] +
                                   log_two_minus_mu - current.log_density) /
                          denominator,
          density_h =
              std::exp(eq.log_hydrogen + current.h[1] - current.log_density) /
              denominator,
          density_he = std::exp(eq.log_helium + log_mu - current.log_density) /
                       denominator;
  const W root_error =
      root_distance + sensitivity_h * saha_error[0] +
      sensitivity_1 * saha_error[1] + sensitivity_2 * saha_error[2] +
      density_h * arithmetic * (1 + std::abs(eq.log_hydrogen)) +
      density_he * arithmetic * (1 + std::abs(eq.log_helium)) +
      arithmetic * (1 + std::abs(ell) + std::abs(current.log_density));
  const std::array<W, 6> relative{
      h1 * (saha_error[0] + root_error),
      h0 * (saha_error[0] + root_error),
      (f1 + f2) * saha_error[1] + f2 * saha_error[2] +
          (f1 + 2 * f2) * root_error,
      f0 * saha_error[1] + f2 * saha_error[2] + std::abs(f2 - f0) * root_error,
      f0 * saha_error[1] + (f0 + f1) * saha_error[2] +
          (2 * f0 + f1) * root_error,
      root_error};
  const std::array<W, 6> logs{current.h[0],  current.h[1],  current.he[0],
                              current.he[1], current.he[2], ell};
  row.admission_status = S::ok;
  auto values = groups(row);
  for (unsigned i = 0; i < values.size(); ++i)
    if (policy.requested_outputs & (1u << i))
      store(*values[i], logs[i], relative[i], i < 5);
}
} // namespace
HydrogenHeliumBatch evaluate_hydrogen_helium_equilibrium(
    std::span<const HydrogenHeliumState> states, HydrogenHeliumPolicy q) {
  HydrogenHeliumBatch out;
  if (!q.requested_outputs || (q.requested_outputs & ~63u) ||
      q.maximum_rows > hard_rows || q.maximum_solves > hard_rows ||
      q.maximum_root_iterations > hard_iterations ||
      q.maximum_charge_evaluations > hard_evaluations ||
      q.maximum_native_bytes > hard_bytes)
    return out;
  if (states.size() > q.maximum_rows || states.size() > hard_rows) {
    out.status = S::work_limit;
    return out;
  }
  std::size_t bytes = sizeof(HydrogenHeliumBatch);
  if (!irred::detail::checked_payload_add(bytes, states.size(),
                                          sizeof(HydrogenHeliumRow)) ||
      bytes > q.maximum_native_bytes) {
    out.status = S::work_limit;
    return out;
  }
  if (std::numeric_limits<W>::digits < 64 ||
      std::numeric_limits<W>::max_exponent < 16384 ||
      std::fegetround() != FE_TONEAREST)
    return out;
  out.rows.resize(states.size());
  out.status = S::ok;
  for (std::size_t i = 0; i < states.size(); ++i) {
    auto &row = out.rows[i];
    row.source = states[i];
    const auto &s = states[i];
    if (!std::isfinite(s.temperature_kelvin) ||
        !std::isfinite(s.hydrogen_nuclei_per_cubic_metre) ||
        !std::isfinite(s.helium_nuclei_per_cubic_metre)) {
      fail(row, S::nonfinite_input, q.requested_outputs);
      continue;
    }
    if (s.temperature_kelvin < 1 || s.temperature_kelvin > 1e5 ||
        !(s.hydrogen_nuclei_per_cubic_metre > 0) ||
        !(s.helium_nuclei_per_cubic_metre > 0)) {
      fail(row, S::outside_domain, q.requested_outputs);
      continue;
    }
    const W kT = boltzmann_constant_joule_per_kelvin * s.temperature_kelvin,
            nq = detail::electron_quantum_density_si(kT),
            maximum_electrons = W(s.hydrogen_nuclei_per_cubic_metre) +
                                2 * W(s.helium_nuclei_per_cubic_metre),
            ratio = maximum_electrons / nq;
    if (ratio * (1 - arithmetic) > .001L) {
      fail(row, S::outside_domain, q.requested_outputs);
      continue;
    }
    if (ratio * (1 + arithmetic) > .001L) {
      fail(row, S::conditioning_budget_exceeded, q.requested_outputs);
      continue;
    }
    if (out.solves == q.maximum_solves) {
      fail(row, S::work_limit, q.requested_outputs);
      continue;
    }
    ++out.solves;
    row.solves = 1;
    const W logq = std::log(nq);
    const std::array<W, 3> binding{
        hydrogen_ionization_energy_ev * electron_volt_joule / kT,
        helium_first_ionization_energy_ev * electron_volt_joule / kT,
        helium_second_ionization_energy_ev * electron_volt_joule / kT};
    const std::array<W, 3> log_saha{
        logq - binding[0], logq + 2 * ln2 - binding[1], logq - binding[2]};
    const std::array<W, 3> log_error{
        arithmetic * (1 + std::abs(logq) + binding[0]),
        arithmetic * (1 + std::abs(logq) + 2 * ln2 + binding[1]),
        arithmetic * (1 + std::abs(logq) + binding[2])};
    const Equation eq{std::log(W(s.hydrogen_nuclei_per_cubic_metre)),
                      std::log(W(s.helium_nuclei_per_cubic_metre)), log_saha};
    solve(row, out, q, eq, std::log(maximum_electrons), log_error);
  }
  return out;
}
} // namespace irred::atomic
