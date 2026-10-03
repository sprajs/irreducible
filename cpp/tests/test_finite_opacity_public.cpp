#include "irred/finite_opacity_source.hpp"
#include "irred/quantities.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

namespace {
using namespace irred::cosmology;
using W=long double;using S=irred::numerics::Status;
void need(bool ok,const char *why){if(!ok){std::cerr<<"FAIL "<<why<<'\n';std::exit(1);}}
// Synthetic axis construction only, sharing the thermal model ancestry. This
// does not act as an independent clock/physics/error reference.
W fixture_age(const ThermalBackground &background,W a) {
  constexpr unsigned n=16384;W sum=0;
  const W c=W(irred::speed_of_light_m_per_s)/1000;
  for(unsigned i=0;i<=n;++i) {
    const auto p=background.scaled_expansion(a*i/n);
    need(p.status==S::ok && p.a4_e2>0,"synthetic fixture clock query");
    sum+=(i==0 || i==n?1:i%2?4:2)*c/background.source().h0_km_s_mpc/std::sqrt(p.a4_e2);
  }
  return a*sum/(3*n);
}
FiniteOpacityRequest request() {
  FiniteOpacityRequest r;r.model={67.4,.0224,.12,2.7255,0,{}};
  r.initial_scale_factor=.001;r.observer_scale_factor=.00101;
  const auto map=map_thermal_physical_model(r.model);need(map.model.has_value(),"synthetic fixture map");
  const auto background=prepare_thermal_background(*map.model);need(background.status()==S::ok,"synthetic fixture background");
  r.opacity.eta_mpc={double(fixture_age(background,r.initial_scale_factor)),
                     double(fixture_age(background,r.observer_scale_factor))};
  r.opacity.differential_opacity_mpc_inverse={.1,.2};
  r.opacity.origin="synthetic-no-electron-history";r.opacity.law_id="emitted-linear-conformal-opacity/v1";
  r.opacity.exact_member_digest="synthetic-fixture-axis-construction-not-content-certification";
  r.k_mpc_inverse={.005,.01};r.source_origin="synthetic-finite-IVP-native-contract";return r;
}
void public_owner() {
  auto r=request();const auto *original=r.opacity.eta_mpc.data();
  FiniteOpacityPolicy trial;
  trial.maximum_attempted_steps=2000;trial.maximum_background_clock_calls=400000;
  trial.maximum_coupled_stage_solves=2000;trial.maximum_destination_writes=100000000;
  trial.maximum_native_bytes=8*1024*1024;
  auto producer=prepare_finite_opacity_source(std::move(r),trial);need(producer.status()==S::ok,"finite owner preparation");
  need(producer.identity() && producer.identity()->original.opacity.eta_mpc.data()==original,"one moved original axis");
  need(producer.identity()->clock.endpoint_clock_consistent && producer.identity()->seeds.size()==2,"same mapped clock and captured seeds");
  const auto *identity=producer.identity();const auto *boundary=&identity->boundary;
  need(boundary->survival_i>0 && boundary->survival_i<1 && boundary->tau_i>0 &&
       boundary->omitted_temperature_absolute_bound[0]>0,"positive surviving boundary and omitted bound");
  auto moved=std::move(producer);need(producer.status()==S::invalid_input && !producer.identity(),"producer move invalidates getters");
  moved=std::move(moved);need(moved.identity()==identity,"self move preserves owner");
  FiniteOpacityPolicy short_work;short_work.maximum_attempted_steps=1;
  const auto refused=moved.produce(short_work);
  need(refused.status==S::work_limit && !refused.source && !refused.source_numerically_admitted &&
       refused.boundary()==boundary && refused.identity.get()==identity &&
       !refused.attempts.empty() && !refused.attempts[0].nodes.empty(),"work refusal preserves original/boundary/prefix");
  need(refused.work.background_clock_calls>=identity->preparation_work.background_clock_calls &&
       refused.work.attempted_steps==1,"original aggregate work is not reset");
  auto result=moved.produce(trial);
  need(result.status==S::conditioning_budget_exceeded && result.source && !result.source_numerically_admitted,
       "complete raw computed source with admission withheld");
  need(result.attempts.size()==9 && result.source->eta_mpc.size()==5 && result.source->t0.size()==10,
       "nine complete trajectories and eta-major original k order");
  for(const auto &attempt:result.attempts) {
    need(attempt.status==S::ok && attempt.nodes.size()==10 && attempt.reached_wavenumbers==2,
         "all hierarchy/time witnesses retained");
    need(attempt.final_temperature_tail.size()==2*(attempt.hierarchy-2) &&
         attempt.final_polarization_tail.size()==2*(attempt.hierarchy-2),"owned full endpoint hierarchy");
  }
  for(const auto &d:result.diagnostics)need(!d.common_background_clock_error && !d.arithmetic_linear_error && !d.source_grid_error,
      "missing integrated errors remain absent");
  need(result.peak_owned_payload_bound && *result.peak_owned_payload_bound<=result.policy.maximum_native_bytes,
       "retained simultaneous payload receipt");
  const auto *channel_buffer=result.source->t0.data();
  auto projection=irred::projection::prepare_continuous_cmb_projection(std::move(*result.source));
  need(projection.status()==S::ok && projection.source()->t0.data()==channel_buffer,
       "one move export of the computed engineering source");
  moved=FiniteOpacitySourceProducer{};
  need(result.identity.get()==identity && result.boundary()==boundary && result.boundary()->survival_i>0,
       "original/boundary lifetime after producer destruction");
  const std::array<unsigned,2> multipoles{2,0};
  const auto projected=irred::projection::project_continuous_cmb(projection,multipoles,
      irred::projection::continuous_temperature);
  need(projected.rows.size()==4,"computed raw source projection rows");
  for(const auto &row:projected.rows)need(!row.e_mode && !row.source_grid_error_estimate && !row.radial_error_estimate,
      "projection masks do not invent producer or radial errors");
}
void invalid_and_profile() {
  FiniteOpacityRequest r;r.model={67.4,.0224,.12,2.7255,0,{}};
  r.initial_scale_factor=.001;r.observer_scale_factor=.002;
  r.opacity.eta_mpc={1,2};r.opacity.differential_opacity_mpc_inverse={0,-1};
  r.opacity.origin="synthetic";r.opacity.law_id="linear";r.opacity.exact_member_digest="fixture";
  r.k_mpc_inverse={.1};r.source_origin="fixture";
  need(prepare_finite_opacity_source(std::move(r)).status()==S::outside_domain,"negative opacity refusal");
  auto duplicate=r;duplicate.opacity.differential_opacity_mpc_inverse={0,0};duplicate.opacity.eta_mpc={1,1};
  need(prepare_finite_opacity_source(std::move(duplicate)).status()==S::outside_domain,"duplicate time refusal");
  auto neutrino=r;neutrino.model.physical_massless_nonphoton_density=.001;
  need(prepare_finite_opacity_source(std::move(neutrino)).status()==S::outside_domain,"different radiation model refused");
  const int saved=std::fegetround();need(std::fesetround(FE_DOWNWARD)==0,"set hostile rounding");
  const auto rejected=prepare_finite_opacity_source(std::move(r));
  need(std::fesetround(saved)==0,"restore rounding");
  need(rejected.status()==S::outside_domain,"distinct strict arithmetic profile refusal");
}
} // namespace
int main(){invalid_and_profile();public_owner();std::cout<<"finite opacity public owner contract passed\n";}
