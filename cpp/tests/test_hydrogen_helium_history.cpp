#include "irred/hydrogen_helium_history.hpp"
#include "../src/hydrogen_helium_rates.hpp"
#include "hydrogen_helium_rate_peer_facts.hpp"
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void need(bool b, const char *message) { ++checks; if (!b) throw std::runtime_error(message); }
double value(const HydrogenHeliumHistoryValue &v) {
  if (v.status != S::ok || !v.value)
    std::cerr << "groupstatus=" << int(v.status) << " error=" << v.absolute_error_estimate << '\n';
  need(v.status == S::ok && v.value.has_value(), "requested group accepted"); return *v.value;
}
HydrogenHeliumHistoryRequest request() {
  return {{67.4, .02237, .12, 2.7255, 1.7e-5, {}}, .19, .015, "synthetic independently supplied nuclei", 2700, 300};
}
} // namespace
int main() {
  try {
    for (const auto &r : hydrogen_helium_rate_peer_facts::rows) {
      const auto h = detail::hhe_rates(r[0], false), he = detail::hhe_rates(r[0], true);
      const auto f = detail::hhe_rhs(.7L, .4L, r[0], 2e9L, 1.5e8L, 1e-13L, 1800);
      const long double actual[]{h.alpha, he.alpha, h.beta, he.beta, h.photo, he.photo, f.value[0], f.value[1]};
      for (unsigned k = 0; k < 8; ++k)
        need(std::abs(actual[k] - r[k + 1]) <= 2e-15L * std::abs(r[k + 1]),
             "independent Decimal110/150 rate and escape facts");
    }
    // Analytic production derivatives face an original centered-difference
    // calculation; the independent history peer uses its own direct equations.
    const long double p = .7L, q = .4L, T = 5000, nH = 2e9L, nHe = 1.5e8L,
        H = 1e-13L, u = 1800, d = 1e-5L;
    auto a = detail::hhe_rhs(p, q, T, nH, nHe, H, u);
    for (unsigned j = 0; j < 2; ++j) {
      auto up = detail::hhe_rhs(p + (j == 0 ? d : 0), q + (j == 1 ? d : 0), T, nH, nHe, H, u);
      auto down = detail::hhe_rhs(p - (j == 0 ? d : 0), q - (j == 1 ? d : 0), T, nH, nHe, H, u);
      for (unsigned i = 0; i < 2; ++i)
        need(std::abs((up.value[i] - down.value[i]) / (2 * d) - a.fraction_derivative[i][j]) <
                 2e-8L * std::abs(a.fraction_derivative[i][j]), "shared charge fraction derivative");
    }
    const auto up = detail::hhe_rhs(p, q, T + .001L, nH, nHe, H, u),
        down = detail::hhe_rhs(p, q, T - .001L, nH, nHe, H, u);
    for (unsigned i = 0; i < 2; ++i)
      need(std::abs((up.value[i] - down.value[i]) / .002L - a.temperature_derivative[i]) <
               2e-8L * std::abs(a.temperature_derivative[i]), "temperature derivative");
    const auto s = request();
    HydrogenHeliumHistoryPolicy fine_policy;
    fine_policy.base_intervals = 16384;
    auto owner = prepare_hydrogen_helium_history(s, fine_policy);
    std::cerr << "prepare=" << int(owner.status()) << " work=" << owner.work().total()
              << " rhs=" << owner.work().rhs_evaluations << '\n';
    need(owner.status() == S::ok && owner.source() && owner.background(), "retained history prepared");
    need(owner.work().total() <= 4000000 && owner.work().initial_charge_evaluations > 0 &&
             owner.work().rhs_evaluations > 0, "charged actual work");
    need(owner.excluded_initial_heiii_activity() && *owner.excluded_initial_heiii_activity() > 0 &&
             owner.maximum_log_excluded_heiii_activity() &&
             *owner.maximum_log_excluded_heiii_activity() <= std::log(1e-12), "explicit omitted-stage witness");
    const double z[]{2700, 2699.99, 2500, 2200, 1900, 1600, 1300, 1100, 800, 600, 300};
    auto batch = owner.evaluate(z, 31);
    need(batch.status == S::ok && batch.rows.size() == 11, "ordered batch");
    for (const auto &r : batch.rows) {
      const double xp = value(r.hydrogen_ionized_fraction), xhe = value(r.helium_singly_ionized_fraction),
          tm = value(r.matter_temperature_kelvin), ne = value(r.electron_number_density_per_cubic_metre),
          op = value(r.thomson_opacity_per_redshift), u = 1 + r.redshift;
      need(xp > 0 && xp < 1 && xhe > 0 && xhe < 1 && tm > 0 && tm <= s.model.tcmb_kelvin * u && op > 0,
           "positive shared-charge physical state");
      need(std::abs(ne - std::pow(u, 3) * (.19 * xp + .015 * xhe)) < 8e-15 * ne,
           "exact declared shared electron closure");
      std::cerr << "z=" << r.redshift << " HII=" << xp << " HeII=" << xhe << " T=" << tm
                << " errors=" << r.hydrogen_ionized_fraction.absolute_error_estimate << ','
                << r.helium_singly_ionized_fraction.absolute_error_estimate << ','
                << r.matter_temperature_kelvin.absolute_error_estimate << '\n';
    }
    need(value(batch.rows.front().matter_temperature_kelvin) == s.model.tcmb_kelvin * 2701,
         "supplied initial temperature");
    // Initial detailed balance is independent of the history stepping.
    const auto &initial = batch.rows.front();
    const long double initialT = value(initial.matter_temperature_kelvin), e = value(initial.electron_number_density_per_cubic_metre);
    const auto h = detail::hhe_rates(initialT, false), he = detail::hhe_rates(initialT, true);
    const long double x = value(initial.hydrogen_ionized_fraction), y = value(initial.helium_singly_ionized_fraction);
    need(std::abs(h.alpha * e * x - h.photo * (1 - x)) < 1e-7L * h.alpha * e * x,
         "initial hydrogen restricted Saha balance");
    need(std::abs(he.alpha * e * y - he.photo * (1 - y)) < 1e-12L * he.alpha * e * y,
         "initial helium factor-four Saha balance");
    for (unsigned mask = 1; mask <= 31; ++mask) {
      auto one = owner.evaluate(std::span(z + 6, 1), mask);
      need(one.status == S::ok, "all nonempty masks");
      const auto &r = one.rows.front();
      const HydrogenHeliumHistoryValue *v[]{&r.hydrogen_ionized_fraction, &r.helium_singly_ionized_fraction,
          &r.electron_number_density_per_cubic_metre, &r.matter_temperature_kelvin, &r.thomson_opacity_per_redshift};
      for (unsigned k = 0; k < 5; ++k) need(bool(v[k]->value) == bool(mask & (1u << k)), "independent output mask");
    }
    const auto before = owner.work().total();
    double invalid[]{1300, 299, std::numeric_limits<double>::quiet_NaN(), 2801, 800};
    auto mixed = owner.evaluate(invalid, 3);
    need(mixed.status == S::ok && mixed.rows.size() == 5 &&
             mixed.rows[1].hydrogen_ionized_fraction.status == S::outside_domain &&
             mixed.rows[2].hydrogen_ionized_fraction.status == S::nonfinite_input &&
             mixed.rows[3].helium_singly_ionized_fraction.status == S::outside_domain &&
             mixed.rows[0].hydrogen_ionized_fraction.value && mixed.rows[4].hydrogen_ionized_fraction.value,
         "invalid queries preserve good rows");
    need(owner.work().total() == before, "queries do no physical work");
    need(owner.evaluate(z, 0).status == S::invalid_input && owner.evaluate(z, 32).status == S::invalid_input,
         "invalid masks");
    need(owner.evaluate(z, 1, 1).status == S::work_limit && owner.evaluate(z, 1, 4096, 1).status == S::work_limit,
         "query resource refusal");
    auto copied = owner;
    owner = HydrogenHeliumHistory{};
    auto moved = std::move(copied);
    need(copied.status() == S::invalid_input && !copied.source() && moved.status() == S::ok,
         "move invalidates original, copy owns state");
    auto *self = &moved; moved = std::move(*self);
    need(moved.status() == S::ok && moved.source()->nuclei_origin == s.nuclei_origin, "self move and source ownership");
    auto wrong = s; wrong.helium_nuclei_today_per_cubic_metre = 0;
    need(prepare_hydrogen_helium_history(wrong).status() == S::outside_domain, "absent species refused");
    wrong = s; wrong.model.species.push_back({.06, 1.95, 2});
    need(prepare_hydrogen_helium_history(wrong).status() == S::outside_domain, "unsupported relic profile explicit");
    wrong = s; wrong.late_redshift = 0;
    need(prepare_hydrogen_helium_history(wrong).status() == S::outside_domain, "late-tail closure absent");
    wrong = s; wrong.initial_redshift = INFINITY;
    need(prepare_hydrogen_helium_history(wrong).status() == S::nonfinite_input, "nonfinite source");
    HydrogenHeliumHistoryPolicy low; low.maximum_total_work = 10;
    auto exhausted = prepare_hydrogen_helium_history(s, low);
    need(exhausted.status() == S::work_limit && exhausted.work().total() <= 10, "original work cap retained");
    low = {}; low.maximum_native_bytes = 1;
    need(prepare_hydrogen_helium_history(s, low).status() == S::work_limit, "preparation byte cap");
    low = {}; low.maximum_fine_intervals = 4;
    need(prepare_hydrogen_helium_history(s, low).status() == S::work_limit, "mesh quota");
    low = {}; low.absolute_fraction_tolerance = -1;
    need(prepare_hydrogen_helium_history(s, low).status() == S::invalid_input, "negative numerical allocation");
    low = {}; low.base_intervals = 1;
    need(prepare_hydrogen_helium_history(s, low).status() == S::invalid_input,
         "coarse curvature needs a resolved stencil");
    auto ordinary = prepare_hydrogen_helium_history(s);
    auto ordinary_row = ordinary.evaluate(std::span(z + 6, 1), 31);
    need(ordinary.status() == S::ok && ordinary_row.status == S::ok &&
             ordinary_row.rows.front().helium_singly_ionized_fraction.status == S::conditioning_budget_exceeded &&
             !ordinary_row.rows.front().helium_singly_ionized_fraction.value &&
             ordinary_row.rows.front().hydrogen_ionized_fraction.value &&
             ordinary_row.rows.front().matter_temperature_kelvin.value,
         "default trace-He diagnostic refusal preserves other groups at unchanged allocation");
    const auto &a8 = ordinary_row.rows.front(), &a16 = batch.rows[6];
    const HydrogenHeliumHistoryValue *c8[]{&a8.hydrogen_ionized_fraction, &a8.electron_number_density_per_cubic_metre,
        &a8.matter_temperature_kelvin, &a8.thomson_opacity_per_redshift};
    const HydrogenHeliumHistoryValue *c16[]{&a16.hydrogen_ionized_fraction, &a16.electron_number_density_per_cubic_metre,
        &a16.matter_temperature_kelvin, &a16.thomson_opacity_per_redshift};
    for (unsigned k = 0; k < 4; ++k)
      need(std::abs(value(*c8[k]) - value(*c16[k])) <= c8[k]->absolute_error_estimate + c16[k]->absolute_error_estimate,
           "explicit native refinement respects retained accepted-group diagnostics");
    low = {}; low.absolute_fraction_tolerance = 1e-30; low.relative_fraction_tolerance = 1e-30;
    auto strict = prepare_hydrogen_helium_history(s, low);
    need(strict.status() == S::ok, "preparation distinct from per-group admission");
    auto partial = strict.evaluate(std::span(z + 6, 1), 9);
    need(partial.status == S::ok &&
             partial.rows.front().hydrogen_ionized_fraction.status == S::conditioning_budget_exceeded &&
             !partial.rows.front().hydrogen_ionized_fraction.value &&
             partial.rows.front().matter_temperature_kelvin.value, "fraction refusal preserves temperature");
    low = {}; low.maximum_native_bytes = *hydrogen_helium_history_payload_bound(
        32768, 0, 0, s.nuclei_origin.size());
    auto capped = prepare_hydrogen_helium_history(s, low);
    need(capped.status() == S::ok && capped.evaluate(std::span(z, 1), 1).status == S::work_limit,
         "query enforces retained byte policy");
    need(!hydrogen_helium_history_payload_bound(SIZE_MAX, 0, 0, 0), "checked payload overflow");
    std::fesetround(FE_UPWARD);
    need(prepare_hydrogen_helium_history(s).status() == S::invalid_input, "unsupported rounding");
    need(moved.evaluate(z, 1).status == S::invalid_input, "prepared query refuses changed rounding");
    std::fesetround(FE_TONEAREST);
    std::cout << "hydrogen_helium_history checks=" << checks << " failures=0\n";
  } catch (const std::exception &e) {
    std::fesetround(FE_TONEAREST); std::cerr << "FAILED " << e.what() << '\n'; return 1;
  }
}
