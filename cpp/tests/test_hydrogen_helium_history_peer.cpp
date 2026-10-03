#include "hydrogen_helium_history_reference.hpp"
#include "irred/hydrogen_helium_history.hpp"
#include <iostream>
#include <stdexcept>
namespace ref = hydrogen_helium_history_reference;
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void need(bool b, const char *message) { ++checks; if (!b) throw std::runtime_error(message); }
double value(const HydrogenHeliumHistoryValue &x) {
  need(x.status == S::ok && x.value.has_value(), "native comparison group accepted"); return *x.value;
}
void compare(ref::W actual, ref::W expected, ref::W allocation, const char *message) {
  const ref::W difference = std::abs(actual - expected);
  if (difference > allocation) std::cerr << message << " delta=" << difference << " allocation=" << allocation << '\n';
  need(difference <= allocation, message);
}
} // namespace
int main() {
  try {
    for (unsigned case_id = 0; case_id < 3; ++case_id) {
      ref::Model model;
      if (case_id == 1) { model.H0 = 60; model.baryon = .015; model.cdm = .08; model.T0 = 2.7;
        model.other = 0; model.nH0 = .1; model.nHe0 = .005; model.initial = 2600; model.late = 600; }
      if (case_id == 2) { model.H0 = 80; model.baryon = .03; model.cdm = .15; model.T0 = 2.75;
        model.other = 3e-5; model.nH0 = .3; model.nHe0 = .035; model.initial = 2800; model.late = 300;
        model.massless_species.push_back({1.95, 2}); }
      HydrogenHeliumHistoryRequest s{{model.H0, model.baryon, model.cdm, model.T0, model.other, {}},
        model.nH0, model.nHe0, "synthetic independent numerical comparison", model.initial, model.late};
      for (const auto &v : model.massless_species) s.model.species.push_back({0, v[0], v[1]});
      std::cerr << "starting_case=" << case_id << " initial=" << model.initial << " late=" << model.late << '\n';
      HydrogenHeliumHistoryPolicy policy; policy.base_intervals = 16384;
      auto native = prepare_hydrogen_helium_history(s, policy);
      need(native.status() == S::ok, "named bounded case prepared");
      const std::vector<ref::W> z{model.initial, model.initial - .01, 2500, 2200, 1900, 1600, 1300, 1100, 800, model.late};
      auto a = ref::integrate(model, z, 1), b = ref::integrate(model, z, 2);
      need(a.rows.size() == z.size() && b.rows.size() == z.size(), "original ordered references complete");
      std::vector<double> queries; for (auto x : z) queries.push_back(double(x));
      auto result = native.evaluate(queries, 31);
      need(result.status == S::ok && result.rows.size() == z.size(), "native ordered queries complete");
      for (std::size_t j = 0; j < z.size(); ++j) {
        const auto &r = result.rows[j];
        const auto &aa = a.rows[j], &bb = b.rows[j];
        need(aa.z == z[j] && bb.z == z[j], "reference source order");
        const ref::W budgets[]{1e-8L + 2e-6L * std::abs(bb.hydrogen),
            1e-8L + 2e-6L * std::abs(bb.helium), 2e-4L + 2e-6L * std::abs(bb.temperature),
            std::pow(1 + z[j], 3) * (model.nH0 + model.nHe0) * 1e-8L + 2e-6L * std::abs(bb.electrons),
            2e-9L + 3e-6L * std::abs(bb.opacity)};
        const ref::W coarse[]{aa.hydrogen, aa.helium, aa.temperature, aa.electrons, aa.opacity},
            refined[]{bb.hydrogen, bb.helium, bb.temperature, bb.electrons, bb.opacity};
        const double observed[]{value(r.hydrogen_ionized_fraction), value(r.helium_singly_ionized_fraction),
            value(r.matter_temperature_kelvin), value(r.electron_number_density_per_cubic_metre), value(r.thomson_opacity_per_redshift)};
        const double diagnostics[]{r.hydrogen_ionized_fraction.absolute_error_estimate,
            r.helium_singly_ionized_fraction.absolute_error_estimate, r.matter_temperature_kelvin.absolute_error_estimate,
            r.electron_number_density_per_cubic_metre.absolute_error_estimate, r.thomson_opacity_per_redshift.absolute_error_estimate};
        for (unsigned k = 0; k < 5; ++k) {
          if (std::abs(observed[k] - refined[k]) > diagnostics[k] + std::abs(coarse[k] - refined[k]) +
              128 * std::numeric_limits<ref::W>::epsilon() * std::abs(refined[k]))
            std::cerr << "diagnostic_case=" << case_id << " z=" << z[j] << " group=" << k
                      << " actual=" << observed[k] << " reference=" << refined[k]
                      << " native_error=" << diagnostics[k] << " reference_refinement="
                      << std::abs(coarse[k] - refined[k]) << '\n';
          compare(coarse[k], refined[k], budgets[k] * .05L, "Radau reference refinement consumes <=5% allocation");
          compare(observed[k], refined[k], budgets[k], "native versus independent stiff history");
          compare(observed[k], refined[k], diagnostics[k] + std::abs(coarse[k] - refined[k]) +
                  128 * std::numeric_limits<ref::W>::epsilon() * std::abs(refined[k]),
                  "retained native diagnostic covers independent discrepancy");
        }
      }
      std::cerr << "case=" << case_id << " work=" << native.work().total() << " reference_steps="
                << a.stats.steps << ',' << b.stats.steps << " reference_rhs=" << a.stats.rhs << ',' << b.stats.rhs << '\n';
    }
    std::cout << "hydrogen_helium_history_peer checks=" << checks << " failures=0\n";
  } catch (const std::exception &e) { std::cerr << "FAILED " << e.what() << '\n'; return 1; }
}
