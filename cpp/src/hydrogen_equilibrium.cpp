#include "irred/hydrogen_equilibrium.hpp"
#include "hydrogen_quantum_density.hpp"
#include "irred/quantities.hpp"
#include "payload_accounting.hpp"
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::atomic {
namespace {
using S = numerics::Status;
using W = long double;
constexpr std::size_t hard_rows = 65536, hard_bytes = 1ull << 30;
constexpr W arithmetic = 64.L * std::numeric_limits<W>::epsilon();
constexpr W relative_admission = 1e-13L;
constexpr W chi = hydrogen_ionization_energy_ev * electron_volt_joule;
void fail(HydrogenRow &r, S s, unsigned mask) noexcept {
  r.admission_status = s;
  if (mask & ionized_fraction)
    r.ionized.status = s;
  if (mask & neutral_fraction)
    r.neutral.status = s;
}
void store(HydrogenFraction &out, W value, W relative) noexcept {
  const double rounded = static_cast<double>(value);
  if (!(value > 0) || !(rounded > 0) || value > 1 || rounded > 1 ||
      !std::isfinite(value) || !std::isfinite(rounded)) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  relative += std::abs(W(rounded) - value) / value + arithmetic;
  if (!std::isfinite(relative) || relative > relative_admission) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out.value = rounded;
  out.relative_arithmetic_estimate = static_cast<double>(relative);
  out.status = S::ok;
}
} // namespace
HydrogenBatch
evaluate_hydrogen_equilibrium(std::span<const HydrogenState> states,
                              HydrogenPolicy q) {
  HydrogenBatch out;
  if (!q.requested_outputs || (q.requested_outputs & ~3u) ||
      q.maximum_rows > hard_rows || q.maximum_solves > hard_rows ||
      q.maximum_native_bytes > hard_bytes)
    return out;
  if (states.size() > hard_rows || states.size() > q.maximum_rows) {
    out.status = S::work_limit;
    return out;
  }
  std::size_t bytes = sizeof(HydrogenBatch);
  if (!irred::detail::checked_payload_add(bytes, states.size(), sizeof(HydrogenRow)) ||
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
        !std::isfinite(s.hydrogen_nuclei_per_cubic_metre)) {
      fail(row, S::nonfinite_input, q.requested_outputs);
      continue;
    }
    if (s.temperature_kelvin < 1 || s.temperature_kelvin > 1e5 ||
        !(s.hydrogen_nuclei_per_cubic_metre > 0)) {
      fail(row, S::outside_domain, q.requested_outputs);
      continue;
    }
    // One owned quantum-density equation serves physical admission and Saha.
    const W kT = boltzmann_constant_joule_per_kelvin * s.temperature_kelvin,
            nq = detail::electron_quantum_density_si(kT),
            density_ratio = W(s.hydrogen_nuclei_per_cubic_metre) / nq;
    if (density_ratio * (1 - arithmetic) > 1e-3L) {
      fail(row, S::outside_domain, q.requested_outputs);
      continue;
    }
    if (density_ratio * (1 + arithmetic) > 1e-3L) {
      fail(row, S::conditioning_budget_exceeded, q.requested_outputs);
      continue;
    }
    if (out.solves == q.maximum_solves) {
      fail(row, S::work_limit, q.requested_outputs);
      continue;
    }
    row.solves = 1;
    ++out.solves;
    const W lnq = std::log(nq),
            lnn = std::log(W(s.hydrogen_nuclei_per_cubic_metre)),
            beta = chi / kT, logk = lnq - beta - lnn;
    // Baseline libm/arithmetic allowance. Sensitivities are d ln x/d ln K
    // = y/(2-x), and |d ln y/d ln K|=x/(2-x), evaluated without 1-x.
    const W log_error = arithmetic * (std::abs(lnq) + beta + std::abs(lnn) + 1);
    W x, y;
    if (logk <= 0) {
      const W rootk = std::exp(logk / 2), root = std::hypot(rootk, W{2});
      x = 2 * rootk / (rootk + root);
      y = 2 / (2 + rootk * rootk + rootk * root);
    } else {
      const W inverse = std::exp(-logk), root = std::sqrt(1 + 4 * inverse);
      x = 2 / (1 + root);
      y = 2 * inverse / (1 + 2 * inverse + root);
    }
    row.admission_status = S::ok;
    if (q.requested_outputs & ionized_fraction)
      store(row.ionized, x, log_error * y / (2 - x));
    if (q.requested_outputs & neutral_fraction)
      store(row.neutral, y, log_error * x / (2 - x));
  }
  return out;
}
} // namespace irred::atomic
