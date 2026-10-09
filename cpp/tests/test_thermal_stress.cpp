#include "irred/thermal_neutrino.hpp"
#include "irred/quantities.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
namespace {
using namespace irred::cosmology;
using S=irred::numerics::Status;
using W=long double;
unsigned checks=0;
void need(bool ok,const char *message) { ++checks; if (!ok) throw std::runtime_error(message); }
void near(W value,W reference,const char *message) {
  need(std::abs(value-reference)<=2e-10L+2e-10L*std::abs(reference),message);
}
W critical() {
  const W pi=std::numbers::pi_v<W>, hbar=6.62607015e-34L/(2*pi),c=299792458.L,ev=1.602176634e-19L;
  const W mpc=648000000000.L*149597870700.L/pi, rate=70*1000/mpc;
  return 3*rate*rate*c*c/(8*pi*6.67430e-11L)/ev*std::pow(hbar*c/ev,3);
}
void analytic() {
  const double factors[]{0,.001,.1,1};
  for (auto model:{ThermalFlatModel{70,1,0,0,0,{}},ThermalFlatModel{70,0,0,1,0,{}},
                   ThermalFlatModel{70,0,0,0,0,{}},ThermalFlatModel{70,.1,.1,.2,.1,{}}}) {
    auto state=prepare_thermal_background(model); need(state.status()==S::ok,"analytic prepare");
    auto result=state.evaluate_stress_energy(factors); need(result.status==S::ok && result.rows.size()==4,"analytic batch");
    const W radiation=model.omega_gamma+W(model.omega_massless_nonphoton),matter=model.omega_b+W(model.omega_cdm),lambda=1-radiation-matter;
    for (const auto &row:result.rows) {
      W a=row.scale_factor,a4=std::pow(a,4),d=radiation+matter*a+lambda*a4,q=radiation/3-lambda*a4;
      if (d==0) { need(row.status==S::conditioning_budget_exceeded && !row.equation_of_state.value,"zero density endpoint refused"); continue; }
      need(row.status==S::ok && row.callbacks==0 && row.species.empty(),"analytic no momentum");
      near(row.scaled_density,d,"analytic scaled density"); near(row.scaled_pressure,q,"analytic scaled pressure");
      need(row.equation_of_state.value && row.deceleration.value,"analytic outputs present");
      near(*row.equation_of_state.value,q/d,"analytic w"); near(*row.deceleration.value,(1+3*q/d)/2,"analytic q");
    }
  }
  auto transition=prepare_thermal_background({70,.3,0,.1,0,{}});
  const double crossing=static_cast<double>(std::pow(.3L/(3*.6L),.25L));
  auto row=transition.evaluate_stress_energy(std::span(&crossing,1)).rows[0];
  near(row.scaled_pressure,0,"pressure cancellation absolute allocation");
  near(*row.equation_of_state.value,0,"zero w absolute allocation");
}
// Independent defining constants and frozen mpmath1.3.0 direct infinite-q
// quadrature60/90digits. No implementation code was copied into the reference.
void massive() {
  auto state=prepare_thermal_background({70,1e-4,2e-4,.05,.25,{{.06,.0002,2}}});
  need(state.status()==S::ok,"massive prepare");
  const double factors[]{0,1./300,.1,1};
  auto result=state.evaluate_stress_energy(factors);
  need(result.status==S::ok && result.rows.size()==4,"massive batch");
  const W factor=2*std::pow(.0002L,4)/(2*std::numbers::pi_v<W>*std::numbers::pi_v<W>*critical());
  const W rho[]{
    5.68219697698347550545901940684113148956744249757331626533562513108163323498158644442377559L,
    6.04564518498587938494993482847105834989123623364985164192184468790340014927826752723655229L,
    54.4781975132089445959504096965791300263163257224343570528725923981828954020959229090659347L,
    540.964487905023651974850009981355610698661719312783465183407491115327904193387256684281758L};
  const W pressure[]{
    1.89406565899449183515300646894704382985581416585777208844520837702721107832719548147459186L,
    1.79179395063083686012993119021955310581949998016371229151753668532579470753205013013092909L,
    0.255012037391274274653278464887215011087519808223257044428722081624381796866411398073149925L,
    0.0259187844188918615737961872805833724916094865066167670288636351765733672683772679568983707L};
  // Filled from immutable reference generation; precision refinement is checked
  // independently before these witnesses are admitted.
  for (unsigned i=0;i<4;++i) {
    auto &row=result.rows[i]; need(row.status==S::ok && row.species.size()==1,"massive per-species stress");
    near(row.species[0].scaled_density/factor,rho[i],"independent FD rho moment");
    near(row.species[0].scaled_pressure/factor,pressure[i],"independent FD pressure moment");
    near(row.scaled_density,state.scaled_expansion(factors[i]).a4_e2,"same retained expansion density");
    need(row.equation_of_state.value && row.deceleration.value,"massive acceleration");
    if (i==3) need(row.callbacks>0,"today pressure requires momentum work");
  }
  const auto cold=evaluate_thermal_moments(1000000);
  need(cold.status==S::ok && cold.rho_moment && cold.pressure_moment,"cold independent moments");
  near(*cold.rho_moment,1803085.35475105686534488082049461230834231039775108693466037L,"cold rho reference");
  near(*cold.pressure_moment,.00000777695816345616318909953190976653704823897066724937418L,"cold pressure absolute reference");
  // Independent derivative diagnostic has own8e-7 discretization allocation.
  // Five-point log-a difference and factor2 refinement, not a2e-10 reference.
  auto derivative=[&](W a,W step) {
    auto logh=[&](W x) { auto p=state.scaled_expansion(a*std::exp(x));
      need(p.status==S::ok,"continuity background"); return std::log(p.a4_e2)/2-2*std::log(a*std::exp(x)); };
    return (-logh(2*step)+8*logh(step)-8*logh(-step)+logh(-2*step))/(12*step);
  };
  for (double a:{.001,.1,.9}) {
    const auto row=state.evaluate_stress_energy(std::span(&a,1)).rows[0];
    const W coarse=-1-derivative(a,.002L),fine=-1-derivative(a,.001L);
    need(std::abs(coarse-fine)<=4e-8L,"continuity refinement<=5percent allocation");
    need(std::abs(*row.deceleration.value-fine)<=8e-7L,"independent log-H continuity diagnostic");
  }
}
void adversarial() {
  auto state=prepare_thermal_background({70,.001,0,.05,.25,{{.06,.0002,2}}});
  const double factors[]{.1,std::numeric_limits<double>::quiet_NaN(),-1,1.1,1};
  auto result=state.evaluate_stress_energy(factors);
  need(result.rows.size()==5 && result.rows[1].status==S::nonfinite_input &&
       result.rows[2].status==S::outside_domain && result.rows[3].status==S::outside_domain &&
       result.rows[4].status==S::ok,"ordered independent failed rows retained");
  auto policy=ThermalStressPolicy{}; policy.thermal.maximum_total_callbacks=1;
  auto refused=state.evaluate_stress_energy(std::span(factors,1),policy);
  need(refused.rows[0].status==S::work_limit && refused.callbacks<=1,"work ceiling failed callbacks count");
  policy={}; policy.absolute_tolerance=1e-30; policy.relative_tolerance=0;
  refused=state.evaluate_stress_energy(std::span(factors,1),policy);
  need(!refused.rows[0].equation_of_state.value,"unattainable budget refused");
  policy={}; policy.thermal.maximum_native_bytes=1;
  need(state.evaluate_stress_energy(factors,policy).status==S::work_limit,"payload admission");
  policy={}; policy.thermal.momentum_method=ThermalMomentumMethod::nested_clenshaw_curtis;
  need(state.evaluate_stress_energy(factors,policy).status==S::invalid_input,"method identity retained");
  auto copy=state; auto moved=std::move(state);
  need(state.status()==S::invalid_input && copy.status()==S::ok && moved.status()==S::ok,"move and copy owners");
  moved=std::move(moved); need(moved.status()==S::ok,"self move");
  std::fesetround(FE_UPWARD); need(copy.evaluate_stress_energy(factors).status==S::invalid_input,"rounding refused"); std::fesetround(FE_TONEAREST);
  need(!thermal_stress_payload_bound(std::numeric_limits<std::size_t>::max(),16),"payload overflow");
}
}
int main() { try { analytic(); massive(); adversarial(); std::cout<<"PASS thermal stress "<<checks<<" checks\n"; }
  catch(const std::exception &e) { std::cerr<<"FAIL "<<e.what()<<'\n'; return 1; } }
