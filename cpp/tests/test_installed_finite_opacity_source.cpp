// Fresh public-header/static-archive consumer. No private equation includes.
#include "irred/finite_opacity_source.hpp"
#include "irred/quantities.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <type_traits>
#include <utility>
namespace {
using namespace irred::cosmology;
using W=long double;using S=irred::numerics::Status;
void need(bool ok,const char *why) {
  if(!ok){std::cerr<<"FAIL installed finite opacity "<<why<<'\n';std::exit(1);}
}
// Synthetic input-axis construction with the SAME thermal ancestry. This is
// installed-engine/lifetime evidence, never an independent physical reference.
W fixture_age(const ThermalBackground &background,W upper) {
  constexpr unsigned panels=16384;W sum=0;
  for(unsigned i=0;i<=panels;++i) {
    const auto p=background.scaled_expansion(upper*i/panels);
    need(p.status==S::ok && p.a4_e2>0,"fixture scaled background");
    const W weight=i==0 || i==panels?1:i%2?4:2;
    sum+=weight*(W(irred::speed_of_light_m_per_s)/1000)/
        background.source().h0_km_s_mpc/std::sqrt(p.a4_e2);
  }
  return upper*sum/(3*panels);
}
FiniteOpacityRequest fixture() {
  FiniteOpacityRequest r;r.model={67.4,.0224,.12,2.7255,0,{}};
  r.initial_scale_factor=.001;r.observer_scale_factor=.00101;
  const auto map=map_thermal_physical_model(r.model);need(map.model.has_value(),"fixture map");
  const auto background=prepare_thermal_background(*map.model);need(background.status()==S::ok,"fixture retained background");
  r.opacity.eta_mpc={double(fixture_age(background,r.initial_scale_factor)),
                     double(fixture_age(background,r.observer_scale_factor))};
  r.opacity.differential_opacity_mpc_inverse={.1,.2};r.k_mpc_inverse={.01};
  r.opacity.origin="installed-synthetic-history";r.opacity.law_id="linear-conformal-opacity/v1";
  r.opacity.exact_member_digest="installed-fixture-construction-not-certification";
  r.source_origin="installed-finite-IVP-engine-contract";return r;
}
} // namespace
int main() {
  static_assert(!std::is_copy_constructible_v<FiniteOpacitySourceProducer>);
  static_assert(!std::is_copy_constructible_v<FiniteOpacitySourceResult>);
  FiniteOpacityRequest invalid;auto refused=prepare_finite_opacity_source(std::move(invalid));
  need(refused.status()!=S::ok && !refused.identity(),"invalid request refused");
  const auto invalid_result=refused.produce();
  need(invalid_result.status!=S::ok && !invalid_result.source && !invalid_result.boundary(),"invalid owner cannot export");
  auto input=fixture();const auto *original_axis=input.opacity.eta_mpc.data();
  FiniteOpacityPolicy policy;policy.maximum_attempted_steps=2000;
  policy.maximum_background_clock_calls=400000;policy.maximum_coupled_stage_solves=2000;
  policy.maximum_destination_writes=100000000;policy.maximum_native_bytes=8*1024*1024;
  auto prepared=prepare_finite_opacity_source(std::move(input),policy);
  need(prepared.status()==S::ok && prepared.identity() &&
       prepared.identity()->original.opacity.eta_mpc.data()==original_axis,"one acquired original axis");
  const auto *identity=prepared.identity();const auto *boundary=&identity->boundary;
  auto moved=std::move(prepared);
  need(prepared.status()==S::invalid_input && !prepared.identity() && moved.identity()==identity,"public move invalidates getters");
  auto result=moved.produce(policy);
  need(result.status==S::conditioning_budget_exceeded && result.source && !result.source_numerically_admitted,
       "complete computed source explicitly unadmitted");
  need(result.identity.get()==identity && result.boundary()==boundary && boundary->survival_i>0 &&
       boundary->survival_i<1 && boundary->omitted_temperature_absolute_bound[0]>0,"same positive boundary owner");
  need(result.attempts.size()==9 && result.source->k_mpc_inverse==identity->original.k_mpc_inverse &&
       result.source->eta_mpc.size()==5 && result.source->eta_mpc.front()==original_axis[0] &&
       result.source->eta_mpc.back()==original_axis[1],"complete original order and finite support");
  std::size_t begun=0;
  for(const auto &attempt:result.attempts) {
    need(attempt.status==S::ok && attempt.nodes.size()==5 && attempt.reached_wavenumbers==1 &&
         attempt.final_core.size()==1 && attempt.denied_step_requests==0,"complete installed trajectory receipts");
    begun+=attempt.attempted_steps;
  }
  need(begun==result.work.attempted_steps && result.work.denied_step_requests==0 &&
       result.peak_owned_payload_bound && *result.peak_owned_payload_bound<=policy.maximum_native_bytes &&
       result.work.background_clock_calls<=policy.maximum_background_clock_calls &&
       result.work.destination_writes<=policy.maximum_destination_writes,"installed work/payload receipt");
  auto short_policy=policy;short_policy.maximum_attempted_steps=1;
  short_policy.maximum_eta_step_mpc=std::min(4.0,(original_axis[1]-original_axis[0])/8);
  const auto stopped=moved.produce(short_policy);
  need(stopped.status==S::work_limit && !stopped.source && stopped.boundary()==boundary &&
       stopped.work.attempted_steps==1 && stopped.work.denied_step_requests==1 &&
       stopped.attempts[0].attempted_steps==1 && stopped.attempts[0].denied_step_requests==1,
       "installed capped refusal retains original identity and causal counters");
  const auto &prefix=stopped.attempts[0];
  need(prefix.nodes.size()==1 && prefix.completed_steps==1 && prefix.reached_wavenumbers==1 &&
       prefix.final_core.size()==1 && prefix.final_eta_mpc.size()==1 && prefix.final_scale_factor.size()==1 &&
       prefix.final_eta_mpc[0]>prefix.nodes[0].eta_mpc &&
       prefix.final_eta_mpc[0]<W(original_axis[0])+(W(original_axis[1])-original_axis[0])/4 &&
       prefix.final_scale_factor[0]>prefix.nodes[0].scale_factor && prefix.final_core[0]!=prefix.nodes[0].core &&
       prefix.final_temperature_tail.size()==prefix.hierarchy-2 &&
       prefix.final_polarization_tail.size()==prefix.hierarchy-2,
       "installed between-node refusal retains matching newer core/time/scale and full tails");
  for(W value:prefix.final_core[0])need(std::isfinite(value),"finite installed refused core");
  for(W value:prefix.final_temperature_tail)need(std::isfinite(value),"finite installed refused temperature tails");
  for(W value:prefix.final_polarization_tail)need(std::isfinite(value),"finite installed refused polarization tails");
  const auto *buffer=result.source->t0.data();
  auto projection=irred::projection::prepare_continuous_cmb_projection(std::move(*result.source));
  need(projection.status()==S::ok && projection.source()->t0.data()==buffer,"one moved raw source buffer");
  moved=FiniteOpacitySourceProducer{};
  need(result.boundary()==boundary && stopped.boundary()==boundary && boundary->survival_i>0,
       "successful/refused receipts outlive producer");
  std::cout<<"PASS installed finite opacity positive computed/unadmitted and refusal contract\n";
}
