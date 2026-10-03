#include "conditional_hydrogen_helium_drag_reference.hpp"
#include "irred/conditional_hydrogen_helium_drag.hpp"
#include <bit>
#include <cfenv>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>

namespace ref = conditional_hydrogen_helium_drag_reference;
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void need(bool condition, const char *message) {
  ++checks;
  if (!condition) throw std::runtime_error(message);
}
void near(ref::W actual, ref::W expected, ref::W allowance, const char *message) {
  const ref::W difference = std::abs(actual - expected);
  if (!(std::isfinite(difference) && difference <= allowance))
    std::cerr << message << " actual=" << actual << " expected=" << expected
              << " difference=" << difference << " allowance=" << allowance << '\n';
  need(std::isfinite(difference) && difference <= allowance, message);
}
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

ConditionalHydrogenHeliumDragRequest request(double late) {
  return {{67.4, .02237, .12, 2.7255, 1.7e-5, {}}, .245,
          1.6735328383153192e-27, 6.646479071583153e-27,
          "NEXT16 v3 exact emitted nuclei/thermal working source",
          "caller-chosen neutral-effective kg masses; no atomic measurement claim",
          2700, late};
}
ConditionalHydrogenHeliumDragPolicy policy() {
  ConditionalHydrogenHeliumDragPolicy p;
  p.history.base_intervals = 16384; p.history.maximum_fine_intervals = 65536;
  p.maximum_total_work = 4000000; p.maximum_native_bytes = 33554432;
  p.history.maximum_total_work = 3999998; p.history.maximum_native_bytes = 33546240;
  p.total_ruler_absolute_tolerance_mpc = 1e-3;
  p.total_ruler_relative_tolerance = 0;
  return p;
}
std::array<ConditionalDragInterval, 3> intervals() {
  return {{{"synthetic-depth-sensitivity", .01, .02,
            "supplied synthetic interval; no tail shape or probability"},
           {"zero-tail-truncated-control", 0, 0,
            "explicit truncated control; no physical negligible-tail assertion"},
           {"collapsed-supplied-depth-control", .01, .01,
            "fixed synthetic supplied-depth control"}}};
}
ref::Source emitted_source(const ConditionalHydrogenHeliumDrag &owner,
                           const ConditionalHydrogenHeliumDragRequest &expected) {
  const auto *original = owner.source(); const auto *snapshot = owner.source_snapshot();
  need(original && snapshot && owner.status() == S::ok, "admitted source snapshot present");
  need(bits(original->model.h0_km_s_mpc) == bits(expected.model.h0_km_s_mpc) &&
       bits(original->model.physical_baryon_density) == bits(expected.model.physical_baryon_density) &&
       bits(original->model.physical_cdm_density) == bits(expected.model.physical_cdm_density) &&
       bits(original->model.tcmb_kelvin) == bits(expected.model.tcmb_kelvin) &&
       bits(original->model.physical_massless_nonphoton_density) ==
           bits(expected.model.physical_massless_nonphoton_density) &&
       bits(original->helium4_neutral_mass_fraction) == bits(expected.helium4_neutral_mass_fraction) &&
       bits(original->hydrogen1_neutral_effective_mass_kg) == bits(expected.hydrogen1_neutral_effective_mass_kg) &&
       bits(original->helium4_neutral_effective_mass_kg) == bits(expected.helium4_neutral_effective_mass_kg) &&
       bits(original->initial_redshift) == bits(expected.initial_redshift) &&
       bits(original->late_redshift) == bits(expected.late_redshift) &&
       original->source_origin == expected.source_origin && original->mass_origin == expected.mass_origin,
       "exact original physical request and provenance retained");
  const auto &H = snapshot->nuclei_today.hydrogen_nuclei,
             &He = snapshot->nuclei_today.helium_nuclei;
  need(H.status == S::ok && He.status == S::ok && H.value && He.value,
       "emitted abundance source values admitted");
  need(bits(*H.value) == snapshot->emitted_nuclei_binary64_bits[0] &&
       bits(*He.value) == snapshot->emitted_nuclei_binary64_bits[1],
       "reference consumes exact emitted nuclei bits");
  const auto &m = snapshot->fixed_mapped_thermal_source;
  need(m.species.empty() && original->model.species.empty(), "reference empty species profile");
  const auto *history = owner.history();
  need(history && history->status() == S::ok && history->source() && history->background(),
       "same accepted history and retained background present");
  const auto &hs = *history->source();
  const auto &background = history->background()->source();
  need(bits(hs.hydrogen_nuclei_today_per_cubic_metre) == bits(*H.value) &&
       bits(hs.helium_nuclei_today_per_cubic_metre) == bits(*He.value) &&
       bits(hs.model.tcmb_kelvin) == bits(original->model.tcmb_kelvin) &&
       bits(hs.initial_redshift) == bits(original->initial_redshift) &&
       bits(hs.late_redshift) == bits(original->late_redshift) &&
       hs.nuclei_origin == original->source_origin,
       "history retains the exact emitted nuclei and prescribed-temperature source");
  need(bits(background.h0_km_s_mpc) == bits(m.h0_km_s_mpc) &&
       bits(background.omega_gamma) == bits(m.omega_gamma) &&
       bits(background.omega_b) == bits(m.omega_b) &&
       bits(background.omega_cdm) == bits(m.omega_cdm) &&
       bits(background.omega_massless_nonphoton) == bits(m.omega_massless_nonphoton) &&
       background.species.empty(), "snapshot is the same fixed retained thermal source");
  need(bits(m.h0_km_s_mpc) == bits(original->model.h0_km_s_mpc), "mapped H0 source identity");
  const double fields[]{m.omega_gamma, m.omega_b, m.omega_cdm, m.omega_massless_nonphoton};
  for (unsigned i = 0; i < 4; ++i)
    need(bits(fields[i]) == bits(snapshot->thermal_source_map[i].emitted_value),
         "ordered scalar witness emitted identity");
  // Copy INPUT values now. Native R0, H, errors, roots and rulers are never used
  // to choose the reference answer, its mesh or any numerical allowance.
  return {m.h0_km_s_mpc, original->model.tcmb_kelvin, m.omega_gamma,
          m.omega_massless_nonphoton, m.omega_b, m.omega_cdm, *H.value, *He.value,
          original->initial_redshift, original->late_redshift};
}
void budget(const ConditionalHydrogenHeliumDrag &owner, const ConditionalDragBatch &batch) {
  const auto *b = owner.budget_diagnostics();
  need(b && b->requested_total_work == 4000000 && b->requested_native_bytes == 33554432 &&
       b->requested_history_work == 3999998 && b->requested_history_native_bytes == 33546240,
       "exact requested policy receipt");
  need(b->history_preparation_attempted && b->served_history_work &&
       b->served_history_native_bytes && b->earned_work_before_history &&
       b->parent_live_bytes_before_history, "attempted history admission metadata complete");
  need(*b->earned_work_before_history <= b->requested_total_work &&
       *b->parent_live_bytes_before_history <= b->requested_native_bytes,
       "parent remaining subtraction admitted");
  need(*b->served_history_work <= b->requested_history_work &&
       *b->served_history_work <= b->requested_total_work - *b->earned_work_before_history &&
       *b->served_history_native_bytes <= b->requested_history_native_bytes &&
       *b->served_history_native_bytes <= b->requested_native_bytes - *b->parent_live_bytes_before_history,
       "served child ceilings respect requested and actual remaining");
  const auto prep = owner.preparation_work().checked_total(), total = batch.work.checked_total();
  need(prep && total && *prep > 0 && *prep <= *total && *total <= 4000000,
       "checked all-attempt batch includes preparation within total");
  const auto p = owner.preparation_work();
  need(p.source_map == 2 && batch.work.source_map == p.source_map &&
       batch.work.background >= p.background && batch.work.momentum >= p.momentum &&
       batch.work.charge >= p.charge && batch.work.rhs >= p.rhs &&
       batch.work.cell >= p.cell && batch.work.prefix >= p.prefix &&
       batch.work.primitive >= p.primitive && batch.work.root >= p.root &&
       batch.work.bound >= p.bound && batch.work.capacity_outer >= p.capacity_outer &&
       batch.work.ruler_outer >= p.ruler_outer,
       "each preparation category survives; exactly two source/map attempts");
  need(batch.retained_payload_bytes <= batch.peak_payload_bytes &&
       batch.peak_payload_bytes <= 33554432, "reported combined payload admission");
}

void mathematical_controls() {
  using ref::W;
  auto weight = [](W z) { return std::pow(1 + z, 6) / 4; };
  const std::array<std::array<W, 2>, 7> cells{{{0, 1}, {0, .25L}, {.25L, .5L}, {.5L, 1},
                                             {0, std::ldexp(1.L, -50)},
                                             {1 - std::ldexp(1.L, -50), 1},
                                             {0, std::ldexp(1.L, -60)}}};
  for (const auto &cell : cells) {
    const W expected = ref::rational_cell_mass(cell[0], cell[1]);
    need(expected > 0, "positive independently factored tiny/full cell");
    near(ref::gauss4(weight, cell[0], cell[1]), expected, 256 * ref::epsilon * expected,
         "independent degree-six Gaussian versus positive exact factorization");
  }
  near(ref::rational_cell_mass(0, 1), 127.L / 28, 128 * ref::epsilon,
       "exact whole polynomial mass");
  near(ref::rational_cell_mass(0, .5L), 2059.L / 3584, 128 * ref::epsilon,
       "exact upper conditional root mass");
  near(ref::rational_cell_mass(0, .25L), 61741.L / 458752, 128 * ref::epsilon,
       "exact lower conditional root mass");
  const W Dlo = 1525.L / 3584, Dhi = 397011.L / 458752;
  near(1 - Dhi, ref::rational_cell_mass(0, .25L), 256 * ref::epsilon,
       "supplied upper depth maps to exact lower polynomial root");
  near(1 - Dlo, ref::rational_cell_mass(0, .5L), 256 * ref::epsilon,
       "supplied lower depth maps to exact upper polynomial root");
  auto primitive = [](W z) { return ref::rational_cell_mass(0, z); };
  const W half = std::ldexp(1.L, -18), EK = half;
  const auto expanded = ref::uniform_bracket(primitive, primitive(.25L), .25L - half,
      .25L + half, .125L, .375L, EK, .25L, .5L);
  near(expanded.root_radius, half, 0, "exact mathematical locator radius retained");
  near(expanded.depth_radius, 4 * half, 0, "same-bracket inverse error is EK/m");
  const W root_ruler = expanded.slope * (expanded.root_radius + expanded.depth_radius),
      history = 256.L / 1048576, ruler_quad = 32.L / 1048576,
      arithmetic = 23.L / 1048576;
  near(root_ruler + history + ruler_quad + arithmetic, 321.L / 1048576, 0,
       "full TOTAL assembly includes locator depth and arithmetic terms");
  need(root_ruler + history + ruler_quad + arithmetic < 1e-3L &&
       root_ruler + history + 1.L / 256 + arithmetic > 1e-3L,
       "same mathematical TOTAL accepts/refuses without weakening allocation");
  bool refused_expansion = false;
  try {
    // Its midpoint shift would fit, but its full locator plus that shift does
    // not. This is the adversarial control for the formerly incomplete check.
    (void)ref::uniform_bracket(primitive, primitive(.25L), .25L - half, .25L + half,
                              .25L - 4 * half, .25L + 4 * half, EK, .25L, .5L);
  } catch (const std::runtime_error &) { refused_expansion = true; }
  need(refused_expansion, "full expanded locator refusal cannot pass midpoint-only containment");
  near(ref::loading_xi(8, 5, 1), .5L, 0, "exact rationalized loading xi control");
  const W tiny_loading_error = std::ldexp(1.L, -std::numeric_limits<W>::digits - 8), A = 1 + 1.5L,
      B = 1 + (1.5L - tiny_loading_error);
  need(std::sqrt(A / B) - 1 == 0 && ref::loading_xi(1.5L, tiny_loading_error, 1) > 0,
       "positive sub-ulp loading correction survives naive subtraction cancellation");
  const W eta = std::ldexp(1.L, -45);
  auto error = [&](W z) {
    const W u = 1 + z, t = z * (1 - z);
    return std::pow(u, 4) * ((2 + eta) * u * t + eta * u * u + t * t) / 4;
  };
  near(ref::gauss5(error, 0, 1), 1037.L / 1260 + 1103 * eta / 224,
       512 * ref::epsilon, "independent degree-eight diagnostic moment");
  need(1 - std::ldexp(1.0, -60) == 1, "trailing rounded public cell is exact zero width");
  const ref::Physics toy({70, 2.7255, 3. / 16, 5. / 16, .25, .25, .2, .02, 1, 0});
  ref::Work work; const ref::Limits limits;
  const W L = ref::ancestor::c / 1000 / 70 * std::sqrt(2.L / 3);
  near(ref::ruler(toy, .5L, 4096, work, limits), L * std::log(5.L / 3), 1e-9L,
       "separate analytic logarithmic ruler upper root");
  near(ref::ruler(toy, .25L, 4096, work, limits), L * std::log(9.L / 5), 1e-9L,
       "separate analytic logarithmic ruler lower root");
  // The polynomial cell and logarithmic ruler are separate mathematical
  // contexts; combining their analytic controls makes no matched HHe source.
}
struct NativeCase { ref::Source source; ConditionalDragBatch batch; };
} // namespace

int main() {
  // Survives every refusal path; no catch or failed gate erases charged work.
  ref::Work work;
  try {
    std::cerr << std::setprecision(21);
    need(std::numeric_limits<ref::W>::digits >= 64 && std::fegetround() == FE_TONEAREST,
         "reference strict wide arithmetic profile");
    std::cerr << "reference_arithmetic_qualification=WITHHELD empirical-work-epsilon-Newton-proxy"
                 " root_complete_reference_gate=WITHHELD\n";
    mathematical_controls();
    const auto tails = intervals(); const auto p = policy();
    std::vector<NativeCase> cases;
    // Finite cases frozen before execution. One matched IVP covers all three
    // late boundaries; no supplied interval becomes a physical tail prediction.
    for (double late : {300., 450., 600.}) {
      const auto input = request(late);
      std::cerr << "native_case late=" << late << " N=16384 fine=65536\n";
      auto owner = prepare_conditional_hydrogen_helium_drag(input, p);
      need(owner.status() == S::ok, "frozen conditional source prepared");
      const auto source = emitted_source(owner, input);
      auto batch = owner.evaluate(tails);
      need(batch.status == S::ok && batch.rows.size() == tails.size(), "ordered native batch complete");
      budget(owner, batch);
      for (std::size_t i = 0; i < tails.size(); ++i) {
        const auto &row = batch.rows[i];
        need(row.source.id == tails[i].id && row.source.origin == tails[i].origin &&
             bits(row.source.lower_depth) == bits(tails[i].lower_depth) &&
             bits(row.source.upper_depth) == bits(tails[i].upper_depth), "exact ordered tail source");
        for (const auto &e : row.endpoints)
          need(e.root_status == S::ok && e.ruler_status == S::ok && e.redshift && e.comoving_ruler_mpc,
               "every requested root/ruler group accepted");
        need(*row.endpoints[0].redshift <= *row.endpoints[1].redshift &&
             *row.endpoints[0].comoving_ruler_mpc >= *row.endpoints[1].comoving_ruler_mpc,
             "conditional inverse and ruler interval orientation");
      }
      if (!cases.empty()) {
        const auto &a = cases.front().source;
        need(bits(source.H0) == bits(a.H0) && bits(source.T0) == bits(a.T0) &&
             bits(source.gamma) == bits(a.gamma) && bits(source.other) == bits(a.other) &&
             bits(source.baryon) == bits(a.baryon) && bits(source.cdm) == bits(a.cdm) &&
             bits(source.nH0) == bits(a.nH0) && bits(source.nHe0) == bits(a.nHe0) &&
             bits(source.initial) == bits(a.initial), "varied boundary retains identical emitted IVP source");
      }
      cases.push_back({source, std::move(batch)});
      if (late == 300) {
        const std::array<ConditionalDragInterval, 2> mixed{{tails[0], {"reversed-control", .03, .02, "synthetic refusal"}}};
        const auto refused = owner.evaluate(mixed);
        need(refused.rows.size() == 2 && refused.rows[1].source.id == mixed[1].id,
             "failed ordered row retained");
        need(refused.rows[0].source.id == tails[0].id &&
             refused.rows[0].source.origin == tails[0].origin &&
             bits(refused.rows[0].source.lower_depth) == bits(tails[0].lower_depth) &&
             bits(refused.rows[0].source.upper_depth) == bits(tails[0].upper_depth),
             "mixed batch preserves accepted source row");
        for (unsigned j = 0; j < 2; ++j) {
          const auto &good = refused.rows[0].endpoints[j],
                     &earlier = cases.back().batch.rows[0].endpoints[j];
          need(good.root_status == S::ok && good.ruler_status == S::ok &&
               good.redshift && good.comoving_ruler_mpc &&
               bits(*good.redshift) == bits(*earlier.redshift) &&
               bits(*good.comoving_ruler_mpc) == bits(*earlier.comoving_ruler_mpc) &&
               bits(good.redshift_numerical_estimate) == bits(earlier.redshift_numerical_estimate) &&
               bits(good.total_ruler_numerical_estimate_mpc) ==
                   bits(earlier.total_ruler_numerical_estimate_mpc),
               "accepted mixed row survives refusal with identical root/ruler diagnostics");
        }
        for (const auto &e : refused.rows[1].endpoints)
          need(e.root_status != S::ok && e.ruler_status != S::ok &&
               !e.redshift && !e.comoving_ruler_mpc, "reversed supplied interval refused without clipping");
        budget(owner, refused);
      }
    }
    using Endpoints = std::vector<std::array<ref::Endpoint, 6>>;
    std::array<Endpoints, 3> reference;
    const ref::Limits limits;
    for (unsigned level = 0; level < 3; ++level) {
      std::cerr << "reference_level=" << level << " refinement=" << (1u << level)
                << " log_step=" << limits.log_step << '\n';
      ref::History history(cases.front().source, 1u << level, work, limits);
      for (const auto &item : cases) {
        std::array<ref::Endpoint, 6> endpoints{};
        for (std::size_t i = 0; i < tails.size(); ++i) {
          std::cerr << "reference_endpoint level=" << level << " late=" << item.source.late
                    << " row=" << tails[i].id << " D_upper=" << tails[i].upper_depth
                    << " D_lower=" << tails[i].lower_depth << '\n';
          endpoints[2 * i] = history.endpoint(tails[i].upper_depth, item.source.late);
          endpoints[2 * i + 1] = history.endpoint(tails[i].lower_depth, item.source.late);
        }
        reference[level].push_back(endpoints);
      }
    }
    for (std::size_t case_id = 0; case_id < cases.size(); ++case_id)
      for (unsigned i = 0; i < 6; ++i) {
        const auto &a = reference[0][case_id][i], &b = reference[1][case_id][i],
                   &fine = reference[2][case_id][i];
        const auto error = ref::complete_error(a, b, fine);
        const auto &native = cases[case_id].batch.rows[i / 2].endpoints[i % 2];
        std::cerr << "endpoint late=" << cases[case_id].source.late << " slot=" << i
                  << " reference_z=" << fine.z << " ruler=" << fine.ruler
                  << " refinement=" << error.history_depth_refinement
                  << " depth=" << error.depth_quadrature << " ruler_quad=" << error.ruler_quadrature
                  << " bracket=" << error.root_bracket << " ALL_arithmetic_proxy=" << error.all_arithmetic
                  << " preliminary_complete=" << error.total()
                  << " native_TOTAL=" << native.total_ruler_numerical_estimate_mpc << '\n';
        need(std::isfinite(error.total()) && error.total() <= ref::endpoint_allowance_mpc,
             "EACH preliminary reference refinement plus arithmetic proxy fits frozen 5e-5 Mpc");
        need(std::abs(b.ruler - fine.ruler) <= std::abs(a.ruler - b.ruler) +
             error.depth_quadrature + error.ruler_quadrature + error.all_arithmetic,
             "reference endpoint refinement resolves rather than hides discrepancy");
        need(std::isfinite(native.total_ruler_numerical_estimate_mpc) &&
             native.total_ruler_numerical_estimate_mpc >= 0 &&
             native.total_ruler_numerical_estimate_mpc <= 1e-3,
             "native unchanged TOTAL admission");
        near(*native.redshift, fine.z, native.redshift_numerical_estimate + error.root_redshift,
             "independent root discrepancy fits complete root estimates");
        near(*native.comoving_ruler_mpc, fine.ruler, 1e-3L,
             "native versus independent ruler frozen allocation");
        near(*native.comoving_ruler_mpc, fine.ruler,
             native.total_ruler_numerical_estimate_mpc + error.total(),
             "native versus independent ruler reported complete estimates");
      }
    std::cerr << "reference steps=" << work.steps << " rhs=" << work.rhs
              << " quadrature=" << work.quadrature << " inverse=" << work.inverse
              << " accepted_relative_corrections=" << work.accepted_relative_corrections << '\n';
    std::cout << "conditional_hydrogen_helium_drag_peer checks=" << checks
              << " failures=0 arithmetic_qualification="
              << (ref::complete_arithmetic_qualified ? "qualified" : "WITHHELD")
              << " root_complete_reference_gate=WITHHELD\n";
  } catch (const std::exception &e) {
    std::cerr << "FAILED " << e.what() << " charged_reference_steps=" << work.steps
              << " rhs=" << work.rhs << " quadrature=" << work.quadrature
              << " inverse=" << work.inverse << " newton=" << work.newton
              << " accepted_relative_corrections=" << work.accepted_relative_corrections << '\n';
    return 1;
  }
}
