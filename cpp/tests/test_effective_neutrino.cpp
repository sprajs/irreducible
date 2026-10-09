#include "irred/effective_neutrino.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace {
using namespace irred::cosmology;
using S=irred::numerics::Status;
unsigned checks=0;
void need(bool value,const char *message) { ++checks; if (!value) throw std::runtime_error(message); }
EffectiveNeutrinoModel anchor() { return {67.4,.02237,.12,2.7255,3.046,{{.06,.71611,1}}}; }
void run() {
  auto source=anchor(); auto state=prepare_effective_neutrino_state(source);
  need(state.status()==S::ok,"source closure preparation");
  auto &d=state.diagnostics();
  need(std::abs(d.remaining_massless_neff-2.0327983725164924L)<2e-14L,"source explicit Neff partition");
  need(std::abs(d.realized_early_neff-3.046)<2e-14,"actual emitted early Neff");
  need(d.realized_early_neff_error_estimate>0 && d.partition_error_estimate>0,"source diagnostics retained");
  need(state.physical_source().species[0].statistical_weight==2,"pair degeneracy to populated states");
  need(d.physical_relic_density_today && *d.physical_relic_density_today>0 && d.omega_lambda,"actual derived present density");
  const double factors[]{0,.001,.1,1}; auto stress=state.evaluate_stress_energy(factors);
  need(stress.status==S::ok && stress.rows.size()==4,"source retained stress batch");
  need(std::abs(*stress.rows[0].equation_of_state.value-1./3)<2e-10,"source analytic radiation limit");
  need(stress.rows[3].callbacks>0,"source present-day pressure not skipped");
  auto request=state.observable_request(1059.95,"synthetic supplied drag","explicit effective thermal realization");
  need(request && request->model.species[0].temperature_today_kelvin==state.physical_source().species[0].temperature_today_kelvin,"consumer same source");
  auto observables=prepare_thermal_observables(*request);
  need(observables.status()==S::ok,"existing consumer accepts new source");
  const double redshift[]{0,.1}; auto distances=observables.evaluate(redshift,early_late_mask(EarlyLateOutput::dh_mpc));
  need(distances.status==S::ok && distances.rows[1].outputs[2].value,"real distance consumer");
  // CLASSv3.3.0 runtime0ceb7a9a; independently executed background owner.
  // Explicit coefficient/phase-space mapping matches native SI/G conventions;
  // qmax80,bins4096/8192 and background10000/20000,tol1e-12 refinement
  // consumes at most1.023percent of the fixed2e-10abs+2e-10rel allocation.
  const double external_a[]{.001,.01,.1,.5,1};
  const long double external[][6]={
    {0.00040560780622385883L, 3.0663011632522549e-05L, 0.075597686144135354L, 0.61339652921620302L, 1.2638388449820502e-05L, 4.1398010850292186e-06L},
    {0.0032331120268209944L, 2.9429024632915586e-05L, 0.0091023832112158871L, 0.51365357481682383L, 1.9537708181255943e-05L, 2.9126644258445462e-06L},
    {0.031630575903582148L, -4.1512064800020891e-05L, -0.0013124030661521933L, 0.49803139540077174L, 0.00014251767432654923L, 4.7498041033398285e-07L},
    {0.2003075209775185L, -0.042792290995300832L, -0.21363297187479854L, 0.17955054218780214L, 0.00070914937118043216L, 9.6070310706134561e-08L},
    {1.0000000000000002L, -0.68507599316579781L, -0.6850759931657977L, -0.52761398974869644L, 0.0014180825451992961L, 4.8052477752429726e-08L}
  };
  const auto matched=state.evaluate_stress_energy(external_a);
  for (unsigned i=0;i<5;++i) {
    const auto &row=matched.rows[i];
    need(row.status==S::ok && row.species.size()==1 && row.equation_of_state.value && row.deceleration.value,"matched CLASS admitted row");
    const long double actual[]{row.scaled_density,row.scaled_pressure,*row.equation_of_state.value,
      *row.deceleration.value,row.species[0].scaled_density,row.species[0].scaled_pressure};
    for (unsigned j=0;j<6;++j)
      need(std::abs(actual[j]-external[i][j])<=2e-10L+2e-10L*std::abs(external[i][j]),"independent matched CLASS stress");
  }
  auto copy=state; source.species[0].mass_ev=1;
  need(copy.source().species[0].mass_ev==.06,"caller storage independent");
  auto moved=std::move(state); need(state.status()==S::invalid_input && moved.status()==S::ok,"source move invalidation");
  moved=std::move(moved); need(moved.status()==S::ok,"source self move");
  auto control=anchor(); control.species.clear(); control.effective_relativistic_species=0;
  auto zero=prepare_effective_neutrino_state(control);
  need(zero.status()==S::ok && zero.physical_source().physical_massless_nonphoton_density==0 &&
       *zero.diagnostics().physical_relic_density_today==0,"zero Neff control");
  control=anchor(); control.species[0].mass_ev=0;
  need(prepare_effective_neutrino_state(control).status()==S::ok,"massless explicit species control");
  control=anchor(); control.species[0].degeneracy=4;
  need(prepare_effective_neutrino_state(control).status()==S::outside_domain,"overpartition no clipping");
  control=anchor(); control.species[0].temperature_ratio=std::numeric_limits<double>::quiet_NaN();
  need(prepare_effective_neutrino_state(control).status()==S::nonfinite_input,"nonfinite species refused");
  control=anchor(); control.physical_cdm_density=2;
  auto excess=prepare_effective_neutrino_state(control);
  need(excess.status()==S::outside_domain && excess.source().physical_cdm_density==2,"failed closure admitted source retained");
  auto policy=ThermalPolicy{}; policy.maximum_total_callbacks=1;
  auto limited=prepare_effective_neutrino_state(anchor(),policy);
  need(limited.status()==S::work_limit && limited.background().preparation_callbacks()<=1 &&
       limited.source().species.size()==1,"failed normalization work/source retained");
  policy={}; policy.maximum_native_bytes=1;
  need(prepare_effective_neutrino_state(anchor(),policy).status()==S::work_limit,"source payload admission");
  policy={}; policy.momentum_method=ThermalMomentumMethod::nested_clenshaw_curtis;
  auto cc=prepare_effective_neutrino_state(anchor(),policy); need(cc.status()==S::ok,"nested CC source preparation");
  ThermalStressPolicy stress_policy; stress_policy.thermal=policy;
  auto cc_result=cc.evaluate_stress_energy(factors,stress_policy);
  need(cc_result.rows.back().deceleration.value && std::abs(*cc_result.rows.back().deceleration.value-
       *stress.rows.back().deceleration.value)<=2e-10,"matched nested method acceleration");
  need(!copy.observable_request(-1,"drag","source") && !copy.observable_request(1,"","source"),"consumer provenance/domain refusal");
  need(!effective_neutrino_payload_bound(std::numeric_limits<std::size_t>::max()),"source payload overflow");
}
}
int main() { try { run(); std::cout<<"PASS effective neutrino "<<checks<<" checks\n"; }
  catch(const std::exception &e) { std::cerr<<"FAIL "<<e.what()<<'\n'; return 1; } }
