#include "hydrogen_helium_supplied_cases.hpp"
#include "hydrogen_helium_supplied_reference.hpp"
#include <bit>
#include <cfenv>
#include <iomanip>
#include <iostream>
#include <string_view>
namespace ref = hydrogen_helium_supplied_reference;
namespace old = hydrogen_helium_history_reference;
namespace cases = hydrogen_helium_supplied_cases;
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0, failures = 0;
void need(bool condition, std::string_view message) {
  ++checks;
  if (!condition) { ++failures; std::cerr << "FAILED " << message << '\n'; }
}
void analytic_temperature_controls() {
  struct ConstantCoupling {
    ref::W g;
    ref::State operator()(ref::W, const ref::State &state) const {
      return {0, 0, -state[2] - g * (state[2] - 1)};
    }
  };
  for (const ref::W g : {0.L, 2.L}) {
    old::Stats stats;
    constexpr ref::W ds = 1e-3L, theta0 = .9L;
    const auto result = old::step(ConstantCoupling{g}, 0, ds, {.7L, .4L, theta0}, stats,
                                 {2e-16L, 4096, true});
    const ref::W equilibrium = g / (1 + g), expected = equilibrium +
        (theta0 - equilibrium) * std::exp(-(1 + g) * ds);
    need(std::abs(result[2] - expected) <= 1e-12L, "independent constant-Compton/adiabatic analytic limit");
    need(result[0] == .7L && result[1] == .4L, "fixed nuclei-fraction analytic control");
  }
}
ref::Model physical(const HydrogenHeliumHistoryRequest &source) {
  ref::Model m;
  m.H0 = source.model.h0_km_s_mpc; m.baryon = source.model.physical_baryon_density;
  m.cdm = source.model.physical_cdm_density; m.T0 = source.model.tcmb_kelvin;
  m.other = source.model.physical_massless_nonphoton_density;
  m.nH0 = source.hydrogen_nuclei_today_per_cubic_metre;
  m.nHe0 = source.helium_nuclei_today_per_cubic_metre;
  m.initial = source.initial_redshift; m.late = source.late_redshift;
  for (const auto &s : source.model.species) m.massless_species.push_back({s.temperature_today_kelvin, s.statistical_weight});
  return m;
}
ref::EmittedSource emitted(const ThermalFlatModel &source) {
  ref::EmittedSource s{source.h0_km_s_mpc, source.omega_gamma, source.omega_massless_nonphoton,
      source.omega_b, source.omega_cdm, {}};
  std::cout << "emitted_background=" << s.H0 << ',' << s.photon << ',' << s.other << ',' << s.baryon << ',' << s.cdm << '\n';
  for (const auto &v : source.species) {
    s.species.push_back({v.mass_ev, v.temperature_today_ev, v.statistical_weight});
    std::cout << "emitted_species=" << v.mass_ev << ',' << v.temperature_today_ev << ',' << v.statistical_weight << '\n';
  }
  return s;
}
ref::Boundary boundary(const HydrogenHeliumSuppliedInitialState &source) {
  return {source.hydrogen_ionized_fraction, source.helium_singly_ionized_fraction, source.matter_temperature_kelvin};
}
bool same_source(const ThermalFlatModel &a, const ref::EmittedSource &b) {
  if (a.h0_km_s_mpc != b.H0 || a.omega_gamma != b.photon || a.omega_massless_nonphoton != b.other ||
      a.omega_b != b.baryon || a.omega_cdm != b.cdm || a.species.size() != b.species.size()) return false;
  for (std::size_t k = 0; k < a.species.size(); ++k)
    if (a.species[k].mass_ev != b.species[k][0] || a.species[k].temperature_today_ev != b.species[k][1] ||
        a.species[k].statistical_weight != b.species[k][2]) return false;
  return true;
}
void source_record(const HydrogenHeliumHistory &owner, const HydrogenHeliumHistoryRequest &source) {
  const auto &m = source.model;
  std::cout << "original_physical_source=" << m.h0_km_s_mpc << ',' << m.physical_baryon_density
            << ',' << m.physical_cdm_density << ',' << m.tcmb_kelvin << ',' << m.physical_massless_nonphoton_density
            << " nuclei_m3=" << source.hydrogen_nuclei_today_per_cubic_metre << ',' << source.helium_nuclei_today_per_cubic_metre
            << " z=" << source.initial_redshift << ',' << source.late_redshift << " origin=" << source.nuclei_origin << '\n';
  for (const auto &item : m.species) std::cout << "original_physical_species=" << item.mass_ev << ','
      << item.temperature_today_kelvin << ',' << item.statistical_weight << '\n';
  if (const auto *maps = owner.thermal_mapping_witnesses()) for (unsigned k = 0; k < 4; ++k) {
    const auto &w = (*maps)[k];
    std::cout << "original_map_witness=" << k << " wide=" << w.wide_value << " operation=" << w.wide_operation_estimate
              << " cast_loss=" << w.measured_absolute_cast_loss << " emitted=" << w.emitted_value << '\n';
  }
  if (owner.excluded_initial_heiii_activity()) std::cout << "native_initial_HeIII_activity=" << *owner.excluded_initial_heiii_activity() << '\n';
  if (owner.maximum_log_excluded_heiii_activity()) std::cout << "native_max_nodal_log_HeIII=" << *owner.maximum_log_excluded_heiii_activity() << '\n';
}
std::array<ref::W, 5> values(const ref::Row &r) {
  return {r.hydrogen, r.helium, r.temperature, r.electrons, r.opacity};
}
std::array<ref::W, 5> allocations(const ref::Row &r, const ref::Model &m) {
  const ref::W u = 1 + r.z;
  return {1e-8L + 2e-6L * std::abs(r.hydrogen), 1e-8L + 2e-6L * std::abs(r.helium),
      2e-4L + 2e-6L * std::abs(r.temperature),
      u * u * u * (ref::W(m.nH0) + m.nHe0) * 1e-8L + 2e-6L * std::abs(r.electrons),
      2e-9L + 3e-6L * std::abs(r.opacity)};
}
struct Refinement {
  std::array<ref::Result, 4> runs;
  std::vector<std::array<ref::W, 5>> errors;
  bool accepted = false;
};
Refinement reference(std::string_view lane, const ref::Model &m, const ref::EmittedSource &source,
                     ref::Boundary input, const std::vector<double> &z,
                     std::optional<ref::State> wide = {}, bool original_input = false) {
  Refinement result;
  const unsigned before = failures;
  for (unsigned r = 0; r < 4; ++r) {
    ref::Policy policy; policy.refinement = r == 0 ? 1 : r == 1 ? 2 : 4;
    if (r == 3) policy.root_tolerance = 2e-18L;
    result.runs[r] = ref::integrate(m, source, input, z, policy, wide, original_input);
    const auto &run = result.runs[r];
    std::cout << "reference_lane=" << lane << " resolution=" << policy.refinement
              << " root=" << policy.root_tolerance << " complete=" << run.complete
              << " rows=" << run.rows.size() << " attempts=" << run.stats.attempts
              << " accepted=" << run.stats.steps << " rejected=" << run.stats.rejected
              << " rhs=" << run.stats.rhs << " newton=" << run.stats.newton
              << " linear_solves=" << run.stats.linear_solves
              << " root_max=" << run.stats.root_correction << " correction_sums="
              << run.stats.accepted_correction_sum[0] << ',' << run.stats.accepted_correction_sum[1]
              << ',' << run.stats.accepted_correction_sum[2] << " failure=" << run.failure << '\n';
    if (run.attempted_temperature)
      std::cout << "reference_attempted_T=" << lane << ',' << r << ','
                << (*run.attempted_temperature)[0] << ',' << (*run.attempted_temperature)[1] << '\n';
    if (run.maximum_log_heiii) std::cout << "reference_attempted_log_HeIII=" << lane << ',' << r << ',' << *run.maximum_log_heiii << '\n';
    need(run.complete && run.rows.size() == z.size(), "complete independent reference including every query");
    need(run.stats.attempts <= policy.maximum_attempts && run.stats.rhs <= policy.maximum_rhs,
         "reference original work caps");
  }
  bool complete = true;
  for (const auto &run : result.runs) complete = complete && run.complete && run.rows.size() == z.size();
  if (!complete) return result; // Failure remains a blocking record, not an accepted subset.
  result.errors.resize(z.size());
  for (std::size_t j = 0; j < z.size(); ++j) {
    const auto a = values(result.runs[0].rows[j]), b = values(result.runs[1].rows[j]),
               c = values(result.runs[2].rows[j]), tight = values(result.runs[3].rows[j]);
    const auto A = allocations(result.runs[2].rows[j], m);
    for (unsigned r = 0; r < 4; ++r) need(result.runs[r].rows[j].z == ref::W(z[j]), "exact restored reference order");
    for (unsigned k = 0; k < 5; ++k) {
      const ref::W d12 = std::abs(a[k] - b[k]), d24 = std::abs(b[k] - c[k]),
          root = std::abs(c[k] - tight[k]),
          arithmetic = 128 * std::numeric_limits<ref::W>::epsilon() * std::abs(c[k]),
          E = std::max(d12, d24) + root + arithmetic;
      result.errors[j][k] = E;
      std::cout << "reference_output=" << lane << ',' << j << ',' << z[j] << ',' << k
                << ',' << a[k] << ',' << b[k] << ',' << c[k] << ',' << tight[k]
                << " differences=" << d12 << ',' << d24 << ',' << root
                << " arithmetic=" << arithmetic << " E=" << E << " allocation=" << A[k] << '\n';
      need(std::isfinite(E) && E <= .05L * A[k], "every reference output consumes <=5% original allocation");
    }
  }
  result.accepted = failures == before;
  return result;
}
void compare_native(std::string_view lane, const HydrogenHeliumHistory &owner,
                    const ref::Model &m, const std::vector<double> &z, const Refinement &reference) {
  const auto batch = owner.evaluate(z, 31);
  std::cout << "native_lane=" << lane << " owner_status=" << int(owner.status()) << " batch_status="
            << int(batch.status) << " rows=" << batch.rows.size() << " work=" << owner.work().total()
            << " import=" << owner.work().initial_boundary_evaluations << " initial_Saha="
            << owner.work().initial_charge_evaluations << " rhs=" << owner.work().rhs_evaluations
            << " model=" << owner.model_identity() << " method=" << owner.method_identity() << '\n';
  need(owner.status() == S::ok && batch.status == S::ok && batch.rows.size() == z.size(), "complete native case");
  if (owner.rate_domain_witness()) {
    const auto &w = *owner.rate_domain_witness();
    if (w.attempted_kelvin_range) std::cout << "native_attempted_T=" << (*w.attempted_kelvin_range)[0] << ',' << (*w.attempted_kelvin_range)[1] << '\n';
    if (w.retained_kelvin_range) std::cout << "native_retained_T=" << (*w.retained_kelvin_range)[0] << ',' << (*w.retained_kelvin_range)[1] << '\n';
    std::cout << "native_invalid_T_attempt=" << w.invalid_temperature_attempted << '\n';
  }
  for (std::size_t j = 0; j < batch.rows.size(); ++j) {
    const auto &row = batch.rows[j];
    // Explicit p,q,T,ne,opacity mapping differs from native mask bit order.
    const HydrogenHeliumHistoryValue *groups[]{&row.hydrogen_ionized_fraction,
        &row.helium_singly_ionized_fraction, &row.matter_temperature_kelvin,
        &row.electron_number_density_per_cubic_metre, &row.thomson_opacity_per_redshift};
    need(row.redshift == z[j], "native original source order");
    for (unsigned k = 0; k < 5; ++k) {
      const auto &g = *groups[k];
      std::cout << "native_output=" << lane << ',' << j << ',' << z[j] << ',' << k
                << " status=" << int(g.status) << " value=";
      if (g.value) std::cout << *g.value; else std::cout << "absent";
      std::cout << " error=" << g.absolute_error_estimate << '\n';
      need(g.status == S::ok && g.value.has_value(), "all five native groups mandatory including initial-layer probes");
      if (!g.value || !reference.accepted) continue;
      const auto expected = values(reference.runs[2].rows[j]);
      const auto A = allocations(reference.runs[2].rows[j], m);
      const ref::W E = reference.errors[j][k], difference = std::abs(ref::W(*g.value) - expected[k]);
      need(std::isfinite(g.absolute_error_estimate) && g.absolute_error_estimate <= A[k], "native own allocation");
      need(difference <= A[k], "native independent discrepancy allocation");
      need(difference <= ref::W(g.absolute_error_estimate) + E, "native diagnostic plus independently admitted reference error");
    }
  }
}
void direct(unsigned id, bool lower, bool upper) {
  auto source = cases::source(id);
  if (lower) { source.initial.matter_temperature_kelvin = 4000; source.initial.source_identity += "#Tm4000"; }
  if (upper) {
    const ref::W Tr = ref::W(source.history.model.tcmb_kelvin) * (1 + ref::W(source.history.initial_redshift));
    source.initial.matter_temperature_kelvin = double(Tr); source.initial.source_identity += "#emittedTr";
    if (ref::W(source.initial.matter_temperature_kelvin) > Tr) {
      const auto refusal = prepare_hydrogen_helium_supplied_history(source);
      need(refusal.status() == S::outside_domain && !refusal.initial_import_witness(), "unclipped upper representation refusal");
      std::cout << "upper_boundary=representation_refused complete_numerical_case=false\n";
      return;
    }
  }
  const auto m = physical(source.history);
  const auto z = cases::queries(id);
  std::cout << "input_identity=" << source.initial.source_identity << " initial=" << source.initial.hydrogen_ionized_fraction
            << ',' << source.initial.helium_singly_ionized_fraction << ',' << source.initial.matter_temperature_kelvin << '\n';
  for (double x : {source.initial.hydrogen_ionized_fraction, source.initial.helium_singly_ionized_fraction, source.initial.matter_temperature_kelvin})
    std::cout << "initial_binary64_bits=" << std::bit_cast<std::uint64_t>(x) << '\n';
  HydrogenHeliumHistoryPolicy policy; policy.base_intervals = 16384;
  auto native = prepare_hydrogen_helium_supplied_history(source, policy);
  source_record(native, source.history);
  need(native.background() != nullptr, "emitted source captured from provenance before reference comparison");
  if (!native.background()) { compare_native("emitted", native, m, z, {}); return; }
  const auto captured = emitted(native.background()->source());
  const auto matched = reference("emitted", m, captured, boundary(source.initial), z);
  compare_native("emitted", native, m, z, matched);
  // Independent original-input construction is informational and has its own
  // complete reference gate. Its difference never enlarges primary allocations.
  const auto original = reference("original-physical-input", m, captured, boundary(source.initial), z, {}, true);
  if (matched.accepted && original.accepted) for (std::size_t j = 0; j < z.size(); ++j) {
    const auto a = values(matched.runs[2].rows[j]), b = values(original.runs[2].rows[j]);
    for (unsigned k = 0; k < 5; ++k) std::cout << "mapping_lane_difference=" << j << ',' << k << ',' << a[k] - b[k] << '\n';
  }
  auto ordinary = prepare_hydrogen_helium_supplied_history(source);
  const auto default_batch = ordinary.evaluate(z, 31);
  std::cout << "default_native_status=" << int(ordinary.status()) << " rows=" << default_batch.rows.size() << '\n';
  // Record every default group; an empirical refusal remains a refusal.
  for (std::size_t j = 0; j < default_batch.rows.size(); ++j) {
    const auto &r = default_batch.rows[j];
    const HydrogenHeliumHistoryValue *g[]{&r.hydrogen_ionized_fraction, &r.helium_singly_ionized_fraction,
        &r.matter_temperature_kelvin, &r.electron_number_density_per_cubic_metre, &r.thomson_opacity_per_redshift};
    for (unsigned k = 0; k < 5; ++k) {
      std::cout << "default_output=" << j << ',' << z[j] << ',' << k << " status=" << int(g[k]->status) << " value=";
      if (g[k]->value) std::cout << *g[k]->value; else std::cout << "absent";
      std::cout << " error=" << g[k]->absolute_error_estimate << '\n';
      if (g[k]->value && matched.accepted) {
        const auto expected = values(matched.runs[2].rows[j]);
        const auto A = allocations(matched.runs[2].rows[j], m);
        const ref::W delta = std::abs(ref::W(*g[k]->value) - expected[k]);
        need(g[k]->status == S::ok && delta <= A[k] &&
             delta <= ref::W(g[k]->absolute_error_estimate) + matched.errors[j][k],
             "every accepted default group faces the unchanged independent gates");
      }
    }
  }
}
void restart(double zi) {
  auto source = cases::source(); source.history.initial_redshift = 2800;
  const auto whole_model = physical(source.history);
  HydrogenHeliumHistoryPolicy policy; policy.base_intervals = 16384;
  auto whole_native = prepare_hydrogen_helium_history(source.history, policy);
  source_record(whole_native, source.history);
  need(whole_native.background() != nullptr, "whole source admitted for restart ancestry");
  if (!whole_native.background()) return;
  const auto captured = emitted(whole_native.background()->source());
  // Scaled-cubic/Saha reference initialization remains independent of native.
  old::Physics initial_physics(whole_model);
  const auto saha = initial_physics.initial();
  const ref::State initial{saha[0], saha[1], ref::W(whole_model.T0) * (1 + ref::W(whole_model.initial))};
  std::vector<double> upstream_z{2800, zi, 2500, 2200, 1900, 1600, 1300, 1100, 800, 600, 300};
  const auto upstream = reference("whole-restricted-Saha", whole_model, captured, {}, upstream_z, initial);
  compare_native("whole-restricted-Saha", whole_native, whole_model, upstream_z, upstream);
  if (!upstream.accepted) { need(false, "restart source unavailable without complete upstream reference"); return; }
  const auto &start = upstream.runs[3].rows[1];
  source.history.initial_redshift = zi;
  source.initial = {double(start.hydrogen), double(start.helium), double(start.temperature),
      "independently refined restricted-Saha upstream explicitly rounded to binary64",
      zi == 2700 ? "NEXT15-controls/v1#restart2700" : "NEXT15-controls/v1#restart2600"};
  const ref::W upstream_values[]{start.hydrogen, start.helium, start.temperature};
  const double rounded_values[]{source.initial.hydrogen_ionized_fraction, source.initial.helium_singly_ionized_fraction,
      source.initial.matter_temperature_kelvin};
  for (unsigned k = 0; k < 3; ++k) std::cout << "restart_boundary=" << k << " wide=" << upstream_values[k]
      << " emitted=" << rounded_values[k] << " cast_loss=" << std::abs(upstream_values[k] - ref::W(rounded_values[k]))
      << " bits=" << std::bit_cast<std::uint64_t>(rounded_values[k]) << '\n';
  const auto model = physical(source.history);
  const std::vector<double> z{zi, zi - .01, 2500, 2200, 1900, 1600, 1300, 1100, 800, 600, 300};
  const auto rounded = reference("restart-exact-emitted", model, captured, boundary(source.initial), z);
  const ref::State wide{start.hydrogen, start.helium, start.temperature};
  const auto unrounded = reference("restart-wide-upstream", model, captured, {}, z, wide);
  auto native = prepare_hydrogen_helium_supplied_history(source, policy);
  source_record(native, source.history);
  need(native.background() && same_source(native.background()->source(), captured),
       "restart retains the same ordered exact emitted thermal coefficients");
  compare_native("restart-exact-emitted", native, model, z, rounded);
  if (rounded.accepted && unrounded.accepted) for (std::size_t j = 0; j < z.size(); ++j) {
    const auto a = values(rounded.runs[2].rows[j]), b = values(unrounded.runs[2].rows[j]);
    for (unsigned k = 0; k < 5; ++k) std::cout << "boundary_input_effect=" << j << ',' << k << ',' << a[k] - b[k]
        << " paired_reference_errors=" << rounded.errors[j][k] << ',' << unrounded.errors[j][k] << '\n';
  }
}
} // namespace
int main(int argc, char **argv) {
  std::cout << std::setprecision(std::numeric_limits<ref::W>::max_digits10);
  std::cout << "reference_method=" << ref::method << " wide_digits=" << std::numeric_limits<ref::W>::digits
            << " wide_max_exponent=" << std::numeric_limits<ref::W>::max_exponent << " rounding=" << std::fegetround()
            << " output_order=p[HII/H],q[HeII/He],Tm[K],ne[m^-3],Thomson_opacity[per_redshift]\n";
  if (argc != 2) { std::cerr << "expected one literal case: B0 B1 B2 T4000 Tr restart2700 restart2600\n"; return 2; }
  try {
    analytic_temperature_controls();
    const std::string_view selected = argv[1];
    if (selected == "B0") direct(0, false, false);
    else if (selected == "B1") direct(1, false, false);
    else if (selected == "B2") direct(2, false, false);
    else if (selected == "T4000") direct(0, true, false);
    else if (selected == "Tr") direct(0, false, true);
    else if (selected == "restart2700") restart(2700);
    else if (selected == "restart2600") restart(2600);
    else { std::cerr << "unknown literal case\n"; return 2; }
    std::cout << "hydrogen_helium_supplied_history_peer case=" << selected << " checks=" << checks
              << " failures=" << failures << " physical_qualification=false\n";
  } catch (const std::exception &error) { std::cerr << "FAILED exception " << error.what() << '\n'; return 1; }
  return failures ? 1 : 0;
}
