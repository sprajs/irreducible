#include "irred/helium_escape_history.hpp"
#include "../src/helium_escape_rates.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>
using namespace irred::cosmology;
using W = long double;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void need(bool b, const char *why) { ++checks; if (!b) throw std::runtime_error(why); }
W stationary(W f, W p, W T, W nH) {
  // Direct independent scalar source normalization; no native quantum helper.
  const W s = 4 * 2.414194e21L * T * std::sqrt(T) / nH * std::exp(-285325.L / T);
  return 2 * s / (p + s + std::sqrt((p + s) * (p + s) + 4 * f * s));
}
double value(const HeliumEscapeHistoryValue &v) {
  need(v.status == S::ok && v.value.has_value(), "requested numerical group accepted");
  return *v.value;
}
}
int main() {
  try {
    std::array<HeliumEscapeDriverKnot, 3> rows{{
        {2700, 1e-13, 3e9, 6000, 1e-5},
        {2400, 1e-13, 3e9, 6000, 1e-5},
        {1900, 1e-13, 3e9, 6000, 1e-5}}};
    const W f = W(.08), p = 1 - W(rows[0].hydrogen_neutral_fraction),
        exact = stationary(f, p, 6000, 3e9);
    std::string identity = "synthetic prescribed constant bath / analytic stationary IVP";
    HeliumEscapeHistoryRequest request{rows, .08, static_cast<double>(exact), identity,
        HeliumEscapeDriverRole::synthetic_prescribed_bath};
    auto history = prepare_helium_escape_history(request);
    std::cout << std::hexfloat << "stationary native_model=" << helium_escape_history_model_id
              << " native_method=" << helium_escape_history_method_id
              << " prepare_status=" << int(history.status()) << " work=" << history.work().total()
              << " emitted_f=" << request.helium_to_hydrogen_nuclei_ratio
              << " emitted_q0=" << request.initial_helium_singly_ionized_fraction
              << " analytic_qeq=" << exact << '\n';
    if (std::numeric_limits<W>::digits < 64 || std::numeric_limits<W>::max_exponent < 16384) {
      need(history.status() == S::invalid_input, "unsupported arithmetic is a refusal");
      std::cout << "helium_escape_history checks=" << checks << " arithmetic-refusal=1\n";
      return 0;
    }
    need(history.status() == S::ok && history.source() &&
        history.witness().source_capture_complete, "owned source/history prepared");
    need(history.work().imported_rows == rows.size() && history.work().total() <= 4000000 &&
        history.work().rate_evaluations == history.witness().attempted_rates,
        "complete actual work and attempted rate accounting");
    need(history.witness().minimum_attempted_temperature_kelvin &&
        history.witness().maximum_attempted_temperature_kelvin &&
        *history.witness().minimum_attempted_temperature_kelvin <= 6000 &&
        *history.witness().maximum_attempted_temperature_kelvin >= 6000 &&
        history.witness().maximum_log_retained_heiii_activity &&
        *history.witness().maximum_log_retained_heiii_activity <= std::log(1e-12),
        "attempted temperature and retained excluded-stage witnesses");
    rows[1].radiation_temperature_kelvin = 100; identity.assign("changed caller storage");
    need(history.source()->knots[1].radiation_temperature_kelvin == 6000 &&
        history.source()->producer_identity.find("analytic stationary") != std::string::npos,
        "immutable captured source survives caller mutation");
    const double z[]{2700, 2650, 2400, 2350, 2200, 1900};
    const W boundary_error = std::abs(W(request.initial_helium_singly_ionized_fraction) - exact),
        exact_ne = 3e9L * p + (3e9L * f) * exact,
        c = irred::speed_of_light_m_per_s,
        sigma = 6.6524587051e-29L;
    const auto batch = history.evaluate(z, 7);
    // Analytic constant-bath flow contracts toward qeq; emitted boundary cast
    // and expected scalar arithmetic are explicit, not claimed exact zeros.
    const W analytic_error = 512 * std::numeric_limits<W>::epsilon() * std::abs(exact),
        qref = boundary_error + analytic_error,
        neref = 3e9L * f * qref + 1024 * std::numeric_limits<W>::epsilon() * exact_ne,
        Aq = 1e-8L + 2e-6L * exact,
        Ane = 3e9L * f * 1e-8L + 2e-6L * exact_ne;
    struct OpacityReference { W value, error, allocation; };
    auto opacity_reference = [&](double redshift) {
      const W expected_op = c * sigma * exact_ne /
          (W(rows[0].hubble_per_second) * (1 + W(redshift))),
          opref = c * sigma * neref /
          (W(rows[0].hubble_per_second) * (1 + W(redshift))) +
          1024 * std::numeric_limits<W>::epsilon() * expected_op;
      return OpacityReference{expected_op, opref, 2e-9L + 3e-6L * expected_op};
    };
    auto record = [&](std::size_t index, const char *group,
                      const HeliumEscapeHistoryValue *v, W reference,
                      W reference_error, W allocation) {
      std::cout << "stationary row=" << index << " requested_z=" << z[index]
                << " emitted_z=";
      if (index < batch.rows.size()) std::cout << batch.rows[index].redshift;
      else std::cout << "absent";
      std::cout << " group=" << group << " status=";
      if (v) std::cout << int(v->status);
      else std::cout << "absent";
      std::cout << " value=";
      if (v && v->value) std::cout << *v->value;
      else std::cout << "absent";
      std::cout << " native_error=";
      if (v) std::cout << v->absolute_error_estimate;
      else std::cout << "absent";
      std::cout << " reference=" << reference << " total_Eref=" << reference_error
                << " allocation=" << allocation << '\n';
    };
    std::cout << "stationary batch_status=" << int(batch.status)
              << " requested_rows=6 actual_rows=" << batch.rows.size()
              << " query_driver_work=" << batch.driver_evaluations << '\n';
    // Emit every original requested group, including missing rows/payloads,
    // before any batch/value/error predicate can throw.
    for (std::size_t i = 0; i < 6; ++i) {
      const auto *row = i < batch.rows.size() ? &batch.rows[i] : nullptr;
      const auto op = opacity_reference(z[i]);
      record(i, "q-per-He", row ? &row->helium_singly_ionized_fraction : nullptr,
             exact, qref, Aq);
      record(i, "ne-per-m3", row ? &row->electron_number_density_per_cubic_metre : nullptr,
             exact_ne, neref, Ane);
      record(i, "qT-per-redshift", row ? &row->thomson_opacity_per_redshift : nullptr,
             op.value, op.error, op.allocation);
    }
    std::cout.flush();
    need(batch.status == S::ok && batch.rows.size() == 6, "ordered coarse query batch");
    for (const auto &row : batch.rows) {
      const W q = value(row.helium_singly_ionized_fraction),
          ne = value(row.electron_number_density_per_cubic_metre),
          op = value(row.thomson_opacity_per_redshift);
      const auto expected = opacity_reference(row.redshift);
      need(qref <= .05L * Aq && neref <= .05L * Ane &&
          expected.error <= .05L * expected.allocation, "complete analytic reference allowance <=5% EACH group");
      need(std::abs(q - exact) <= row.helium_singly_ionized_fraction.absolute_error_estimate +
          qref, "analytic stationary HISTORY control");
      need(std::abs(ne - exact_ne) <= row.electron_number_density_per_cubic_metre.absolute_error_estimate +
          neref, "analytic shared electron history");
      need(std::abs(op - expected.value) <= row.thomson_opacity_per_redshift.absolute_error_estimate +
          expected.error, "analytic opacity history");
      need(q > 0 && q < 1 && ne > 0 && op > 0, "positive perHe/charge/opacity outputs");
    }
    for (unsigned mask = 1; mask <= 7; ++mask) {
      auto one = history.evaluate(std::span(z + 2, 1), mask);
      need(one.status == S::ok && one.rows.size() == 1, "all valid masks");
      const auto &r = one.rows[0];
      const HeliumEscapeHistoryValue *v[]{&r.helium_singly_ionized_fraction,
          &r.electron_number_density_per_cubic_metre, &r.thomson_opacity_per_redshift};
      for (unsigned k = 0; k < 3; ++k)
        need(bool(v[k]->value) == bool(mask & (1u << k)), "independent output mask semantics");
    }
    const double invalid[]{2400, 1899, std::numeric_limits<double>::quiet_NaN(), 2801, 2200};
    auto mixed = history.evaluate(invalid, 7);
    need(mixed.status == S::ok && mixed.rows.size() == 5 &&
        mixed.rows[0].helium_singly_ionized_fraction.value &&
        mixed.rows[1].helium_singly_ionized_fraction.status == S::outside_domain &&
        mixed.rows[2].thomson_opacity_per_redshift.status == S::nonfinite_input &&
        mixed.rows[3].electron_number_density_per_cubic_metre.status == S::outside_domain &&
        mixed.rows[4].thomson_opacity_per_redshift.value, "individual query refusals preserve valid rows");
    need(history.evaluate(z, 0).status == S::invalid_input &&
        history.evaluate(z, 8).status == S::invalid_input &&
        history.evaluate(z, 7, 1).status == S::work_limit &&
        history.evaluate(z, 7, 4096, 1).status == S::work_limit, "mask/shape/query-payload refusal");
    need(!helium_escape_history_payload_bound(SIZE_MAX, 2, 1, 0) &&
        !helium_escape_history_payload_bound(4, SIZE_MAX, 1, 0), "payload products checked");
    auto copy = history;
    auto moved = std::move(copy);
    need(copy.status() == S::invalid_input && !copy.source() &&
        !copy.witness().source_capture_complete && moved.status() == S::ok, "copy/move owns immutable source");
    auto *self = &moved; moved = std::move(*self);
    need(moved.source() && moved.evaluate(z, 7).status == S::ok, "self move preserves history");
    HeliumEscapeHistory assigned; assigned = moved;
    need(assigned.status() == S::ok && assigned.source()->producer_identity == moved.source()->producer_identity,
        "copy replacement ownership");
    rows[1].radiation_temperature_kelvin = 6000;
    request.producer_identity = "synthetic refusal source";
    auto low = HeliumEscapeHistoryPolicy{}; low.maximum_total_work = rows.size();
    auto exhausted = prepare_helium_escape_history(request, low);
    need(exhausted.status() == S::work_limit && exhausted.source() &&
        exhausted.witness().source_capture_complete && exhausted.work().total() == rows.size() &&
        exhausted.work().driver_evaluations == 0 && exhausted.work().rate_evaluations == 0,
        "remaining work checked after complete import before numerical setup");
    low = {}; low.maximum_total_work = rows.size() + 4 * low.base_intervals + 2;
    auto partial = prepare_helium_escape_history(request, low);
    std::cout << "partial-work status=" << int(partial.status()) << " source_complete="
              << partial.witness().source_capture_complete << " imported=" << partial.work().imported_rows
              << " driver=" << partial.work().driver_evaluations << " rates=" << partial.work().rate_evaluations
              << " attempted_rates=" << partial.witness().attempted_rates << '\n';
    need(partial.status() == S::work_limit && partial.source() &&
        partial.witness().source_capture_complete && partial.work().rate_evaluations > 0 &&
        partial.witness().minimum_attempted_temperature_kelvin && partial.evaluate(z, 7).rows.empty(),
        "source/work/attempted witness survive partial numerical failure");
    low = {}; low.maximum_native_bytes = 1;
    need(prepare_helium_escape_history(request, low).status() == S::work_limit, "preparation payload refusal");
    low = {}; low.maximum_fine_intervals = 4;
    need(prepare_helium_escape_history(request, low).status() == S::work_limit, "mesh profile refusal");
    low = {}; low.absolute_fraction_tolerance = -1;
    need(prepare_helium_escape_history(request, low).status() == S::invalid_input, "invalid allocation refused");
    low = {}; low.absolute_fraction_tolerance = 1e-30; low.relative_fraction_tolerance = 1e-30;
    auto strict = prepare_helium_escape_history(request, low);
    auto independent = strict.evaluate(std::span(z + 2, 1), 5);
    need(strict.status() == S::ok && independent.status == S::ok &&
        !independent.rows[0].helium_singly_ionized_fraction.value &&
        independent.rows[0].helium_singly_ionized_fraction.status == S::conditioning_budget_exceeded &&
        independent.rows[0].thomson_opacity_per_redshift.value,
        "q diagnostic refusal leaves raw inherited opacity diagnostic independently gated");
    auto wrong = request; wrong.initial_helium_singly_ionized_fraction = .5;
    rows[0].radiation_temperature_kelvin = 8000;
    auto domain = prepare_helium_escape_history(wrong);
    need(domain.status() == S::outside_domain && domain.source() &&
        domain.witness().source_capture_complete, "captured physical input refusal");
    rows[0].radiation_temperature_kelvin = 6000;
    wrong = request; wrong.initial_helium_singly_ionized_fraction = 1e-8;
    auto trace = prepare_helium_escape_history(wrong);
    need(trace.status() == S::outside_domain && trace.witness().attempted_rates > 0 &&
        trace.witness().minimum_attempted_temperature_kelvin, "no hidden trace-floor freeze or dropped interval");
    rows[1].redshift = rows[0].redshift;
    need(prepare_helium_escape_history(request).status() == S::invalid_input, "source order preserved/refused");
    rows[1].redshift = 2400;
    wrong = request; wrong.helium_to_hydrogen_nuclei_ratio = INFINITY;
    need(prepare_helium_escape_history(wrong).status() == S::nonfinite_input, "nonfinite source refused");
    std::fesetround(FE_UPWARD);
    need(prepare_helium_escape_history(request).status() == S::invalid_input &&
        history.evaluate(z, 7).status == S::invalid_input, "arithmetic contract at preparation and query");
    std::fesetround(FE_TONEAREST);
    // Analytic derivative faces centered finite differences, not its own formula.
    const W q = .4L, d = 1e-6L;
    const auto r = detail::helium_escape_rate(1e-13L, 3e9L, 6000, 1e-5L, .08L, q),
        up = detail::helium_escape_rate(1e-13L, 3e9L, 6000, 1e-5L, .08L, q + d),
        down = detail::helium_escape_rate(1e-13L, 3e9L, 6000, 1e-5L, .08L, q - d);
    const W finite_difference = (up.value - down.value) / (2 * d),
        derivative_error = std::abs(finite_difference - r.fraction_derivative),
        derivative_scale = std::max(std::abs(r.fraction_derivative), r.absolute_rate_scale),
        derivative_allocation = 2e-8L * derivative_scale;
    std::cout << "derivative q=" << q << " step=" << d << " center_rate=" << r.value
              << " upper_rate=" << up.value << " lower_rate=" << down.value
              << " analytic=" << r.fraction_derivative << " finite_difference=" << finite_difference
              << " absolute_difference=" << derivative_error << " absolute_scale=" << derivative_scale
              << " allocation=" << derivative_allocation << '\n';
    std::cout.flush();
    need(derivative_error <= derivative_allocation, "analytic q derivative");
    std::cout << "helium_escape_history checks=" << checks << " failures=0\n";
  } catch (const std::exception &e) {
    std::fesetround(FE_TONEAREST); std::cerr << "FAILED " << e.what() << '\n'; return 1;
  }
}
