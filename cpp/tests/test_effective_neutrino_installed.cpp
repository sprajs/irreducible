// Separately compiled against the installed archive and public headers too.
#include <irred/effective_neutrino.hpp>
#include <iostream>
int main() {
  using namespace irred::cosmology;
  auto state=prepare_effective_neutrino_state({67.4,.02237,.12,2.7255,3.046,{{.06,.71611,1}}});
  const double a[]{0,.1,1}; const auto result=state.evaluate_stress_energy(a);
  if (state.status()!=irred::numerics::Status::ok || result.rows.size()!=3 ||
      !result.rows[2].deceleration.value || result.rows[2].callbacks==0 ||
      !state.diagnostics().physical_relic_density_today) return 1;
  const auto request=state.observable_request(1059.95,"synthetic supplied drag","effective FD realization");
  if (!request) return 1;
  const auto observations=prepare_thermal_observables(*request);
  const double z[]{0,.1}; const auto distances=observations.evaluate(z,early_late_mask(EarlyLateOutput::dh_mpc));
  if (!distances.rows[1].outputs[2].value) return 1;
  std::cout<<"PASS installed effective neutrino consumer\n";
}
