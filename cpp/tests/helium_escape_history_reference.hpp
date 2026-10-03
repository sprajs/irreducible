#pragma once
#include "irred/helium_escape_history.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
// Original two-stage, order3 Radau-IIA reference. It consumes finite source data
// only; no native rates, quantum/charge/interpolation/root helpers or trajectory.
// Physical coefficients and system libm ancestry are intentionally shared.
namespace helium_escape_reference {
using W = long double;
using Knot = irred::cosmology::HeliumEscapeDriverKnot;
struct State { W H, nH, T, xHI; };
struct Rate {
  W tau, inverse, width, tc, enhancement, escape, down, saha, ne, value, scale;
};
inline Rate direct_rate(State c, W f, W q) {
  const W s0 = 4 * 2.414194e21L * c.T * std::sqrt(c.T) / c.nH,
      s = s0 * std::exp(-285325.L / c.T),
      ys = std::exp(46090.L / c.T) / s0,
      yp = 3 * std::exp(39101.L / c.T) / s0,
      inverse = 9.15776e28L * c.H / c.nH / c.xHI,
      width = 1.976e6L / (1 - std::exp(-6989.L / c.T)) +
          6.03e6L / (std::exp(19754.L / c.T) - 1) +
          1.06e8L / (std::exp(21539.L / c.T) - 1) +
          2.18e6L / (std::exp(28496.L / c.T) - 1) +
          3.37e7L / (std::exp(29224.L / c.T) - 1) +
          1.04e6L / (std::exp(32414.L / c.T) - 1) +
          1.51e7L / (std::exp(32781.L / c.T) - 1),
      tau = 4.277e-14L * c.nH / c.H * f * (1 - q),
      pi = std::numbers::pi_v<W>, tc = width * tau / (4 * pi * pi * inverse),
      enhancement = std::sqrt(1 + pi * pi * tc) + 7.74L * tc / (1 + 70 * tc),
      B = 1 - std::exp(-1.023e-7L * tau),
      escape = enhancement / tau + B *
          (.964525L * std::exp(2947.L / c.T) -
           enhancement * std::exp(-6.14e13L / inverse)) / tau,
      down = 50.94L * ys + 1.7989e9L * yp * escape,
      xe = 1 - c.xHI + f * q,
      value = down * ((1 - q) * s - q * xe) / c.H,
      scale = down * ((1 - q) * s + q * xe) / c.H;
  return {tau, inverse, width, tc, enhancement, escape, down, s, c.nH * xe, value, scale};
}
inline State driver(std::span<const Knot> rows, W ell) {
  const W z = std::expm1(-ell);
  auto low = std::lower_bound(rows.begin(), rows.end(), z,
      [](const Knot &k, W x) { return W(k.redshift) > x; });
  std::size_t i = low - rows.begin();
  if (!i) i = 1;
  if (i == rows.size()) i = rows.size() - 1;
  const auto &a = rows[i - 1], &b = rows[i];
  const W ea = -std::log1p(W(a.redshift)), eb = -std::log1p(W(b.redshift)),
      t = (ell - ea) / (eb - ea);
  auto geometric = [&](double x, double y) {
    if (x == y) return W(x);
    return W(x) * std::exp(t * std::log(W(y) / W(x)));
  };
  return {geometric(a.hubble_per_second, b.hubble_per_second),
          geometric(a.hydrogen_nuclei_per_cubic_metre, b.hydrogen_nuclei_per_cubic_metre),
          geometric(a.radiation_temperature_kelvin, b.radiation_temperature_kelvin),
          geometric(a.hydrogen_neutral_fraction, b.hydrogen_neutral_fraction)};
}
// Passive last-tested evidence, in per-He q units and dimensionless ln(a).
// Only a fully evaluated full/two-half rejection can populate this fixed record.
struct RejectedTrial {
  std::size_t query_index;
  double query_redshift;
  unsigned retry_index;
  W maximum_step, source_span;
  std::size_t call_cap;
  W old_ell, target_ell, next_split_ell, h, old_q, old_accumulated;
  W q_full, ell_after_full, accumulated_after_full,
      full_root_correction, full_arithmetic_charge;
  W q_half1, ell_after_half1, accumulated_after_half1,
      half1_root_correction, half1_arithmetic_charge;
  W q_half2, ell_after_half2, accumulated_after_half2,
      half2_root_correction, half2_arithmetic_charge;
  W refinement_component, charged_full_difference, charged_half_difference,
      complete_local, local_allowance, epsilon;
  std::size_t calls_before_trial, calls_after_full, calls_after_half1, calls_after_half2;
};
static_assert(sizeof(RejectedTrial) <= 1024);
static_assert(sizeof(std::optional<RejectedTrial>) <= 1024);
struct StepCharge { W root_correction, arithmetic_charge; };
struct Result {
  bool complete = false;
  std::size_t calls = 0;
  std::vector<W> q, numerical_arithmetic_error;
  std::string refusal;
  std::optional<RejectedTrial> rejected_trial;
};
// Adaptive step-doubled Radau with a fixed maximum step, independently refined;
// split every emitted source knot
// and output coordinate. No native grid or extrapolated states are read.
inline Result integrate(std::span<const Knot> rows, W f, W q0,
                         std::span<const double> z, W maximum_step,
                         std::size_t call_cap) {
  Result out;
  if (rows.size() < 2 || !(maximum_step > 0) || !call_cap) return out;
  out.q.reserve(z.size()); out.numerical_arithmetic_error.reserve(z.size());
  W ell = -std::log1p(W(rows.front().redshift)), q = q0, accumulated = 0;
  const W eps = std::numeric_limits<W>::epsilon(),
      source_span = std::log1p(W(rows.front().redshift)) - std::log1p(W(rows.back().redshift));
  auto rate = [&](W t, W x) {
    if (out.calls >= call_cap) throw std::runtime_error("independent reference call limit");
    ++out.calls;
    const auto c = driver(rows, t);
    const auto r = direct_rate(c, f, x);
    if (!(x > 0 && x < 1) || !(r.tau >= 100) || !(f * x >= 1e-6L) ||
        !std::isfinite(r.value)) throw std::runtime_error("independent reference profile refusal");
    return r;
  };
  constexpr W A[2][2]{{5.L / 12, -1.L / 12}, {3.L / 4, 1.L / 4}},
      c[2]{1.L / 3, 1};
  auto step = [&](W h) {
    W x[2]{q, q};
    for (unsigned iteration = 0; iteration < 32; ++iteration) {
      Rate r[2]; W derivative[2];
      for (unsigned j = 0; j < 2; ++j) {
        const W t = ell + c[j] * h, delta = 1e-5L * std::max(x[j], W(1e-4));
        r[j] = rate(t, x[j]);
        const auto plus = rate(t, x[j] + delta), minus = rate(t, x[j] - delta);
        derivative[j] = (plus.value - minus.value) / (2 * delta);
      }
      W R[2], J[2][2];
      for (unsigned i = 0; i < 2; ++i) {
        R[i] = x[i] - q - h * (A[i][0] * r[0].value + A[i][1] * r[1].value);
        for (unsigned j = 0; j < 2; ++j) J[i][j] = (i == j ? 1 : 0) - h * A[i][j] * derivative[j];
      }
      const W det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
      if (!std::isfinite(det) || !(std::abs(det) > 64 * eps))
        throw std::runtime_error("independent reference unresolved Radau stages");
      const W d[2]{(-J[1][1] * R[0] + J[0][1] * R[1]) / det,
                   (J[1][0] * R[0] - J[0][0] * R[1]) / det},
          correction = std::max(std::abs(d[0]), std::abs(d[1]));
      if (correction <= 64 * eps * std::max(x[0], x[1])) {
        // The endpoint is stage2 (stiffly accurate Radau), with accumulated
        // root/finite arithmetic diagnostics. This is not a rigorous certificate.
        const W arithmetic_charge = 1024 * eps *
            (std::abs(q) + std::abs(x[1]) + h * (r[0].scale + r[1].scale));
        accumulated += correction + arithmetic_charge;
        q = x[1]; ell += h; return StepCharge{correction, arithmetic_charge};
      }
      W damping = 1; bool accepted = false;
      for (unsigned k = 0; k < 64; ++k) {
        if (x[0] + damping * d[0] > 0 && x[0] + damping * d[0] < 1 &&
            x[1] + damping * d[1] > 0 && x[1] + damping * d[1] < 1) {
          x[0] += damping * d[0]; x[1] += damping * d[1]; accepted = true; break;
        }
        damping /= 2;
      }
      if (!accepted) throw std::runtime_error("independent reference root refusal");
    }
    throw std::runtime_error("independent reference iteration limit");
  };
  bool final_refinement_exhaustion = false;
  try {
    for (double query : z) {
      const W target = -std::log1p(W(query));
      if (!std::isfinite(query) || query < rows.back().redshift || query > rows.front().redshift ||
          target < ell) return out;
      while (ell < target) {
        W next = target;
        for (const auto &row : rows) {
          const W boundary = -std::log1p(W(row.redshift));
          if (boundary > ell && boundary < next) next = boundary;
        }
        W h = std::min(maximum_step, next - ell);
        bool accepted = false;
        for (unsigned retry = 0; retry < 64; ++retry) {
          const W old_ell = ell, old_q = q, old_error = accumulated;
          const std::size_t calls_before_trial = out.calls;
          const auto full_charge = step(h);
          const W full_q = q, full_error = accumulated, full_ell = ell;
          const std::size_t calls_after_full = out.calls;
          ell = old_ell; q = old_q; accumulated = old_error;
          const auto half1_charge = step(h / 2);
          const W half1_q = q, half1_error = accumulated, half1_ell = ell;
          const std::size_t calls_after_half1 = out.calls;
          const auto half2_charge = step(h / 2);
          const W refinement_component = std::abs(q - full_q) / 7,
              charged_half_difference = accumulated - old_error,
              charged_full_difference = full_error - old_error,
              complete_local = refinement_component + charged_half_difference + charged_full_difference,
              allowance = (1e-12L + 1e-10L * std::max(old_q, q)) * h / source_span;
          if (complete_local <= allowance) {
            out.rejected_trial.reset();
            accumulated = old_error + complete_local;
            // Coordinate endpoint arithmetic stays in the reference error.
            const W endpoint = old_ell + h;
            if (!(endpoint > old_ell)) throw std::runtime_error("independent reference coordinate refusal");
            ell = endpoint; accepted = true; break;
          }
          // Capture the actually tested h and decision operands before rollback.
          // These copies add no rate/driver evaluations or acceptance predicate.
          out.rejected_trial = RejectedTrial{
              out.q.size(), query, retry, maximum_step, source_span, call_cap,
              old_ell, target, next, h, old_q, old_error,
              full_q, full_ell, full_error, full_charge.root_correction, full_charge.arithmetic_charge,
              half1_q, half1_ell, half1_error, half1_charge.root_correction, half1_charge.arithmetic_charge,
              q, ell, accumulated, half2_charge.root_correction, half2_charge.arithmetic_charge,
              refinement_component, charged_full_difference, charged_half_difference,
              complete_local, allowance, eps,
              calls_before_trial, calls_after_full, calls_after_half1, out.calls};
          ell = old_ell; q = old_q; accumulated = old_error; h /= 2;
        }
        if (!accepted) {
          final_refinement_exhaustion = true;
          throw std::runtime_error("independent reference refinement refusal");
        }
        if (next - ell <= 16 * eps * std::max(W(1), std::abs(next))) {
          accumulated += 1024 * eps * std::abs(q); ell = next;
        }
      }
      out.q.push_back(q); out.numerical_arithmetic_error.push_back(accumulated);
    }
    out.rejected_trial.reset();
    out.complete = true;
  } catch (const std::runtime_error &e) {
    out.refusal = e.what();
    if (!final_refinement_exhaustion) out.rejected_trial.reset();
    // Keep original partial outputs/counter as failure evidence; caller must
    // require complete before any reference or output acceptance.
  }
  return out;
}
} // namespace helium_escape_reference
