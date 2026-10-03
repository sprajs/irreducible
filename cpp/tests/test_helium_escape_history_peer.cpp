#include "irred/helium_escape_history.hpp"
#include "../src/helium_escape_rates.hpp"
#include "helium_escape_history_reference.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using W = long double;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void need(bool b, const char *why) { ++checks; if (!b) throw std::runtime_error(why); }
void compare_scalar(W actual, W reference, W scale, W reference_error) {
  const W allocation = 2e-13L * scale;
  std::cout << std::hexfloat << "scalar-comparison actual=" << actual << " reference=" << reference
            << " scale=" << scale << " reference_error=" << reference_error
            << " allocation=" << allocation << '\n';
  need(reference_error <= .05L * allocation, "complete scalar reference numerical allowance <=5%");
  need(std::abs(actual - reference) <= allocation + reference_error, "independent scalar algebra");
}
void emit_reference_failure_prefix(unsigned fixture, const char *level,
                                   std::span<const double> queries,
                                   const helium_escape_reference::Result &result) {
  std::cerr << std::hexfloat << "reference-failure-prefix fixture=" << fixture << " level=" << level
            << " complete=" << result.complete << " refusal=" << result.refusal
            << " calls=" << result.calls << " witness_present=" << result.rejected_trial.has_value()
            << " query_size=" << queries.size() << " q_size=" << result.q.size()
            << " numerical_arithmetic_error_size=" << result.numerical_arithmetic_error.size()
            << " numerical_arithmetic_error_meaning=accumulated-local-refinement-root-arithmetic-diagnostic\n";
  for (std::size_t i = 0; i < std::max(result.q.size(), result.numerical_arithmetic_error.size()); ++i) {
    std::cerr << "reference-prefix-row fixture=" << fixture << " level=" << level << " query_index=" << i
              << " query_redshift=";
    if (i < queries.size()) std::cerr << queries[i]; else std::cerr << "absent";
    std::cerr << " q=";
    if (i < result.q.size()) std::cerr << result.q[i]; else std::cerr << "absent";
    std::cerr << " numerical_arithmetic_error=";
    if (i < result.numerical_arithmetic_error.size()) std::cerr << result.numerical_arithmetic_error[i];
    else std::cerr << "absent";
    std::cerr << " accepted_complete_reference=false\n";
  }
  if (!result.rejected_trial) return;
  const auto &r = *result.rejected_trial;
  std::cerr << "reference-rejected-trial fixture=" << fixture << " level=" << level
            << " query_index=" << r.query_index << " query_redshift=" << r.query_redshift
            << " retry_index=" << r.retry_index << " maximum_step=" << r.maximum_step
            << " source_span=" << r.source_span << " call_cap=" << r.call_cap << '\n'
            << "reference-rejected-input fixture=" << fixture << " level=" << level
            << " old_ell=" << r.old_ell << " target_ell=" << r.target_ell
            << " next_split_ell=" << r.next_split_ell << " h=" << r.h
            << " old_q=" << r.old_q << " old_accumulated=" << r.old_accumulated << '\n'
            << "reference-rejected-full fixture=" << fixture << " level=" << level
            << " q_full=" << r.q_full << " ell_after_full=" << r.ell_after_full
            << " accumulated_after_full=" << r.accumulated_after_full
            << " full_root_correction=" << r.full_root_correction
            << " full_arithmetic_charge=" << r.full_arithmetic_charge << '\n'
            << "reference-rejected-half1 fixture=" << fixture << " level=" << level
            << " q_half1=" << r.q_half1 << " ell_after_half1=" << r.ell_after_half1
            << " accumulated_after_half1=" << r.accumulated_after_half1
            << " half1_root_correction=" << r.half1_root_correction
            << " half1_arithmetic_charge=" << r.half1_arithmetic_charge << '\n'
            << "reference-rejected-half2 fixture=" << fixture << " level=" << level
            << " q_half2=" << r.q_half2 << " ell_after_half2=" << r.ell_after_half2
            << " accumulated_after_half2=" << r.accumulated_after_half2
            << " half2_root_correction=" << r.half2_root_correction
            << " half2_arithmetic_charge=" << r.half2_arithmetic_charge << '\n'
            << "reference-rejected-decision fixture=" << fixture << " level=" << level
            << " refinement_component=" << r.refinement_component
            << " charged_full_difference=" << r.charged_full_difference
            << " charged_half_difference=" << r.charged_half_difference
            << " complete_local=" << r.complete_local << " local_allowance=" << r.local_allowance
            << " epsilon=" << r.epsilon << '\n'
            << "reference-rejected-calls fixture=" << fixture << " level=" << level
            << " calls_before_trial=" << r.calls_before_trial << " calls_after_full=" << r.calls_after_full
            << " calls_after_half1=" << r.calls_after_half1 << " calls_after_half2=" << r.calls_after_half2 << '\n';
}
}
int main() {
  try {
    if (std::numeric_limits<W>::digits < 64 || std::numeric_limits<W>::max_exponent < 16384) {
      std::cout << "helium_escape_history_peer arithmetic unsupported; no numerical acceptance\n"; return 0;
    }
    std::size_t scalar_reference_rates = 0;
    for (W T : {4300.L, 6000.L, 7800.L}) {
      std::cout << std::hexfloat << "scalar-source T=" << T
                << " H=" << 1e-13L << " nH=" << 3e9L << " xHI=" << 1e-5L
                << " f=" << .08L << " q=" << .4L << '\n';
      const auto n = detail::helium_escape_rate(1e-13L, 3e9L, T, 1e-5L, .08L, .4L);
      const auto r = helium_escape_reference::direct_rate({1e-13L, 3e9L, T, 1e-5L}, .08L, .4L);
      ++scalar_reference_rates;
      const W a[]{n.tau, n.continuum_inverse_seconds, n.width_per_second, n.continuum_depth,
          n.enhancement, n.effective_escape, n.downward_per_second, n.saha, n.electron_density},
          b[]{r.tau, r.inverse, r.width, r.tc, r.enhancement, r.escape, r.down, r.saha, r.ne};
      for (unsigned i = 0; i < 9; ++i)
        compare_scalar(a[i], b[i], std::abs(b[i]), 1024 * std::numeric_limits<W>::epsilon() * std::abs(b[i]));
      // Signed net rates use an absolute positive reaction scale, including at zero.
      compare_scalar(n.value, r.value, r.scale, 1024 * std::numeric_limits<W>::epsilon() * r.scale);
      const W nH_cgs = 3e9L * 1e-6L,
          I_cgs = 9.15776e22L * 1e-13L / nH_cgs / 1e-5L,
          tau_cgs = 4.277e-8L * nH_cgs / 1e-13L * (.08L - .08L * .4L),
          per_H = r.down * ((.08L - .08L * .4L) * r.saha -
                           (.08L * .4L) * (1 - 1e-5L + .08L * .4L)) / 1e-13L;
      compare_scalar(n.continuum_inverse_seconds, I_cgs, I_cgs,
                     1024 * std::numeric_limits<W>::epsilon() * I_cgs);
      compare_scalar(n.tau, tau_cgs, tau_cgs,
                     1024 * std::numeric_limits<W>::epsilon() * tau_cgs);
      compare_scalar(.08L * n.value, per_H, .08L * r.scale,
                     1024 * std::numeric_limits<W>::epsilon() * .08L * r.scale);
      // Proper-time/ln(a)/redshift and solver-minus-z roles remain distinct.
      const W dq_dz = -n.value / 2701, wrapper_dq_dz = -per_H / 2701 / .08L;
      compare_scalar(dq_dz, wrapper_dq_dz, r.scale / 2701,
                     2048 * std::numeric_limits<W>::epsilon() * r.scale / 2701);
    }
    {
      const W T = 6000, f = .08L, p = 1 - 1e-5L,
          s = 4 * 2.414194e21L * T * std::sqrt(T) / 3e9L * std::exp(-285325.L / T),
          qeq = 2 * s / (p + s + std::sqrt((p + s) * (p + s) + 4 * f * s));
      const auto equilibrium = detail::helium_escape_rate(1e-13L, 3e9L, T, 1e-5L, f, qeq),
          below = detail::helium_escape_rate(1e-13L, 3e9L, T, 1e-5L, f, qeq / 2),
          above = detail::helium_escape_rate(1e-13L, 3e9L, T, 1e-5L, f, 3 * qeq / 2);
      compare_scalar(equilibrium.value, 0, equilibrium.absolute_rate_scale,
                     1024 * std::numeric_limits<W>::epsilon() * equilibrium.absolute_rate_scale);
      need(below.value > 0 && above.value < 0 && equilibrium.fraction_derivative < 0,
           "detailed-balance stationary root and attracting derivative");
      const auto free = detail::helium_escape_rate(1e-13L, 3e9L, T, 1e-120L, f, .4L);
      const W B = -std::expm1(-1.023e-7L * free.tau),
          limit = (std::exp(-1.023e-7L * free.tau) +
                   B * .964525L * std::exp(2947.L / T)) / free.tau;
      compare_scalar(free.effective_escape, limit, limit,
                     1024 * std::numeric_limits<W>::epsilon() * limit);
      const W q = 1 - 101 * 1e-12L / (4.277e-14L * 5e8L * .04L);
      const auto thick = detail::helium_escape_rate(1e-12L, 5e8L, T, 1e-5L, .04L, q);
      const W x = 1.023e-7L * thick.tau, exact_B = -std::expm1(-x),
          series_B = x - x * x / 2 + x * x * x / 6,
          series_error = std::exp(std::abs(x)) * x * x * x * x / 24 +
                         1024 * std::numeric_limits<W>::epsilon() * exact_B;
      compare_scalar(exact_B, series_B, exact_B, series_error);
      need(thick.tau >= 100 && thick.tau < 102 && .04L * q >= 1e-6L,
           "thin intercombination series retains allowed-line thick regime");
    }
    std::size_t total_reference_calls = 0;
    for (unsigned fixture = 0; fixture < 2; ++fixture) {
      std::array<HeliumEscapeDriverKnot, 3> rows{{
          {2700, 1e-13, 3e9, fixture ? 6200. : 6000., 1e-5},
          {2400, 1e-13, 3e9, 6000, 1e-5},
          {1900, 1e-13, 3e9, fixture ? 5700. : 6000., 1e-5}}};
      HeliumEscapeHistoryRequest request{rows, .08, .02,
          fixture ? "synthetic geometric temperature / emitted conditional IVP" :
                    "synthetic constant bath / nonstationary emitted conditional IVP",
          HeliumEscapeDriverRole::synthetic_prescribed_bath};
      auto owner = prepare_helium_escape_history(request);
      std::cout << "fixture=" << fixture << " native_model=" << helium_escape_history_model_id
                << " native_method=" << helium_escape_history_method_id
                << " reference_method=original-adaptive-Radau3-step-doubling/v1"
                << " prepare_status=" << int(owner.status()) << " work=" << owner.work().total() << '\n';
      for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto &r = rows[i];
        std::cout << std::hexfloat << "driver fixture=" << fixture << " row=" << i
                  << " z=" << r.redshift << " H=" << r.hubble_per_second
                  << " nH=" << r.hydrogen_nuclei_per_cubic_metre
                  << " Tr=" << r.radiation_temperature_kelvin << " xHI=" << r.hydrogen_neutral_fraction
                  << " f=" << request.helium_to_hydrogen_nuclei_ratio
                  << " q0=" << request.initial_helium_singly_ionized_fraction << '\n';
      }
      need(owner.status() == S::ok, "native synthetic trajectory prepared");
      // Every source knot, geometric cell midpoint and both endpoints appear.
      const double z[]{2700, std::sqrt(2701. * 2401.) - 1, 2400,
                              std::sqrt(2401. * 1901.) - 1, 1900};
      const W span = std::log(2701.L / 1901.L);
      auto coarse = helium_escape_reference::integrate(rows, W(.08), W(.02), z, span / 8192, 998976);
      auto fine = helium_escape_reference::integrate(rows, W(.08), W(.02), z, span / 16384, 998976);
      total_reference_calls += coarse.calls + fine.calls;
      std::cout << "reference fixture=" << fixture << " coarse_complete=" << coarse.complete
                << " fine_complete=" << fine.complete << " coarse_calls=" << coarse.calls
                << " fine_calls=" << fine.calls << '\n';
      if (!coarse.complete || !fine.complete) {
        std::cerr << "reference refusal coarse=" << coarse.refusal << " fine=" << fine.refusal
                  << " calls=" << coarse.calls << "," << fine.calls << '\n';
        emit_reference_failure_prefix(fixture, "coarse", z, coarse);
        emit_reference_failure_prefix(fixture, "fine", z, fine);
      }
      need(coarse.complete && fine.complete && fine.q.size() == 5 && coarse.q.size() == 5,
           "complete independently refined Radau trajectory without dropped states");
      const auto native = owner.evaluate(z, 7);
      need(native.status == S::ok && native.rows.size() == 5, "complete coarse native output batch");
      for (std::size_t i = 0; i < 5; ++i) {
        const W ell = -std::log1p(W(z[i]));
        const auto state = helium_escape_reference::driver(rows, ell);
        const W q = fine.q[i], qref = std::abs(q - coarse.q[i]) +
                fine.numerical_arithmetic_error[i] + coarse.numerical_arithmetic_error[i],
            ne = state.nH * (1 - state.xHI + W(.08) * q),
            ne_error = state.nH * W(.08) * qref +
                2048 * std::numeric_limits<W>::epsilon() * (ne + state.nH),
            coefficient = 299792458.L * 6.6524587051e-29L / (state.H * (1 + W(z[i]))),
            opacity = coefficient * ne,
            opacity_error = coefficient * ne_error + 2048 * std::numeric_limits<W>::epsilon() * opacity,
            reference[]{q, ne, opacity}, errors[]{qref, ne_error, opacity_error},
            allocations[]{1e-8L + 2e-6L * std::abs(q),
                state.nH * W(.08) * 1e-8L + 2e-6L * std::abs(ne),
                2e-9L + 3e-6L * std::abs(opacity)};
        const auto &row = native.rows[i];
        const HeliumEscapeHistoryValue *values[]{&row.helium_singly_ionized_fraction,
            &row.electron_number_density_per_cubic_metre, &row.thomson_opacity_per_redshift};
        for (unsigned k = 0; k < 3; ++k) {
          constexpr const char *ids[]{"q-per-He", "ne-per-m3", "qT-per-redshift"};
          std::cout << std::hexfloat << "history-comparison fixture=" << fixture << " row=" << i
                    << " output=" << ids[k] << " z=" << row.redshift
                    << " native_status=" << int(values[k]->status) << " native=";
          if (values[k]->value) std::cout << *values[k]->value; else std::cout << "absent";
          std::cout << " reference=" << reference[k] << " native_error=" << values[k]->absolute_error_estimate
                    << " complete_reference_error=" << errors[k] << " allocation=" << allocations[k]
                    << " coarse_q=" << coarse.q[i] << " fine_q=" << fine.q[i]
                    << " coarse_numerical_error=" << coarse.numerical_arithmetic_error[i]
                    << " fine_numerical_error=" << fine.numerical_arithmetic_error[i] << '\n';
          need(errors[k] <= .05L * allocations[k], "COMPLETE per-output reference numerical error <=5%");
          need(values[k]->status == S::ok && values[k]->value,
               "every separately requested native group accepted");
          need(std::abs(W(*values[k]->value) - reference[k]) <=
              W(values[k]->absolute_error_estimate) + errors[k], "independent history discrepancy fits diagnostics");
        }
      }
    }
    need(total_reference_calls + scalar_reference_rates <= 4000000,
         "whole peer independent rate-request bound, including scalar references");
    std::cout << "helium_escape_history_peer checks=" << checks
              << " reference_calls=" << total_reference_calls << " failures=0\n";
  } catch (const std::exception &e) { std::cerr << "FAILED " << e.what() << '\n'; return 1; }
}
