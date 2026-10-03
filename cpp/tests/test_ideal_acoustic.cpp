#include "irred/ideal_acoustic.hpp"
#include "../src/ideal_acoustic_equations.hpp"
#include "../src/ideal_acoustic_response.hpp"
#include "../src/ideal_acoustic_transport.hpp"
#include "ideal_acoustic_transport_controls.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>
#include <utility>
using namespace irred::cosmology;
using S=irred::numerics::Status;
using W=long double;
void need(bool value,const char *message) {
  if (!value) { std::cerr<<"FAIL "<<message<<'\n'; std::exit(1); }
}
void close(W actual,W expected,W error,const char *message) {
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      !std::isfinite(error) || error<0 || std::abs(actual-expected)>error) {
    std::cerr<<static_cast<double>(actual)<<" vs "<<static_cast<double>(expected)<<'\n';
    need(false,message);
  }
}
void print_failure(const IdealAcousticAttempt &trial) {
  if (!trial.failure) { std::cout<<" failure=null"; return; }
  const auto &f=*trial.failure;
  std::cout<<" failure_stage="<<static_cast<unsigned>(f.stage)
           <<" rk_stage="<<f.rk_stage<<" committed_a=";
  if (f.last_committed_scale_factor) std::cout<<*f.last_committed_scale_factor;
  else std::cout<<"null";
  std::cout<<" coefficient_a=";
  if (f.attempted_coefficient_scale_factor) std::cout<<*f.attempted_coefficient_scale_factor;
  else std::cout<<"null";
  std::cout<<" radius_assembly=";
  if (!f.radius_assembly) { std::cout<<"null"; return; }
  const auto &r=*f.radius_assembly;
  std::cout<<"coordinate:"<<static_cast<unsigned>(r.coordinate)
           <<",channel:"<<static_cast<unsigned>(r.channel)
           <<",arithmetic_status:"<<static_cast<int>(r.arithmetic_status);
  const auto scalar=[](const char *id,const auto &value) {
    std::cout<<','<<id<<':';
    if (value) std::cout<<*value; else std::cout<<"null";
  };
  scalar("provisional",r.provisional_radius);
  scalar("assembly_allowance",r.local_assembly_allowance);
  scalar("completed",r.completed_radius);
  std::cout<<",finite:"<<r.completed_finite<<",normal_or_zero:"<<r.completed_normal_or_zero
           <<",nonnegative:"<<r.completed_nonnegative;
}
void failure_witness_limits(const IdealAcousticTransfer &owner) {
  const std::array<double,1> k{1e-4},a{.01};
  IdealAcousticPolicy age;
  age.maximum_age_evaluations=0;
  const auto age_failure=owner.evaluate(k,a,1,age);
  need(age_failure.status==S::work_limit && age_failure.trajectories.size()==1 &&
       age_failure.trajectories[0].attempts_started==1,"early age refusal shape");
  const auto &age_trial=age_failure.trajectories[0].attempts[0];
  need(age_trial.failure && age_trial.failure->stage==IdealAcousticFailureStage::initial_age &&
       age_trial.failure->rk_stage==0 && !age_trial.failure->last_committed_scale_factor &&
       !age_trial.failure->attempted_coefficient_scale_factor &&
       !age_trial.failure->radius_assembly &&
       age_failure.evaluation_work.background_evaluations==0 &&
       age_failure.evaluation_work.rhs_evaluations==0,"unavailable initial epoch stays absent without callbacks");
  IdealAcousticPolicy storage;
  storage.maximum_state_element_writes=0;
  const auto storage_failure=owner.evaluate(k,a,1,storage);
  need(storage_failure.status==S::work_limit && storage_failure.trajectories.size()==1 &&
       storage_failure.trajectories[0].attempts_started==1,"initial storage refusal shape");
  const auto &storage_trial=storage_failure.trajectories[0].attempts[0];
  need(storage_trial.failure &&
       storage_trial.failure->stage==IdealAcousticFailureStage::initial_state_storage &&
       !storage_trial.failure->last_committed_scale_factor &&
       !storage_trial.failure->radius_assembly &&
       storage_trial.failure->attempted_coefficient_scale_factor==W(owner.source()->initial_scale_factor) &&
       storage_trial.work.state_element_writes==0 && storage_trial.work.rhs_evaluations==0,
       "storage refusal retains attempted coefficient epoch before a state exists");
  IdealAcousticPolicy rhs;
  rhs.maximum_rhs_per_wavenumber=0;
  const auto rhs_failure=owner.evaluate(k,a,1,rhs);
  need(rhs_failure.status==S::work_limit && rhs_failure.trajectories.size()==1 &&
       rhs_failure.trajectories[0].attempts_started==1,"initial RHS refusal shape");
  const auto &rhs_trial=rhs_failure.trajectories[0].attempts[0];
  need(rhs_trial.failure && rhs_trial.failure->stage==IdealAcousticFailureStage::rk_nominal_rhs &&
       rhs_trial.failure->rk_stage==1 &&
       !rhs_trial.failure->radius_assembly &&
       rhs_trial.failure->last_committed_scale_factor==W(owner.source()->initial_scale_factor) &&
       rhs_trial.failure->attempted_coefficient_scale_factor==rhs_trial.failure->last_committed_scale_factor &&
       rhs_trial.work.rhs_evaluations==0 && !rhs_failure.rows[0].outputs[0].computed,
       "first RK refusal retains matching committed and attempted epochs without derivative execution");
}
IdealAcousticRequest fixture() {
  return {{70,.02,.10,2.7,0,{}},1e-10,"synthetic:ideal-acoustic-original150249"};
}
void equation_limits() {
  // Formal radiation limit: all retained regular-series terms must satisfy
  // their independent leading derivatives. This is an equation control, not
  // an admitted baryon-free photon source.
  const W q=1e-8L,p0=-2.L/3;
  std::array<W,5> y{q/4,q*q/144,q/36,p0*(.5L-q/120),p0*(1-q/30)},dy{};
  const irred::cosmology::detail::IdealAcousticCoefficients radiation{
      q,-1,4.L/3,4.L/3,0,1.L/3,0};
  irred::cosmology::detail::ideal_acoustic_derivative(y,radiation,dy);
  close(dy[0],q/2,1e-17L,"radiation Delta leading coefficient");
  close(dy[1],q*q/36,1e-30L,"entropy x4 derivative");
  close(dy[2],q/18,1e-17L,"radiation relative velocity leading coefficient");
  close(dy[3],-p0*q/60,1e-17L,"radiation CDM velocity wave correction");
  close(dy[4],-p0*q/15,1e-17L,"radiation potential wave correction");
  // Exact formal EdS growing solution, Delta=(2/5)zeta*x², phi=-3/5.
  const W x2=.1L;
  y={.4L*x2,0,0,-.4L,-.6L};
  const irred::cosmology::detail::IdealAcousticCoefficients dust{x2,-.5L,1,0,0,0,0};
  irred::cosmology::detail::ideal_acoustic_derivative(y,dust,dy);
  close(dy[0],y[0],1e-18L,"EdS comoving growth");
  close(dy[3],0,1e-18L,"EdS velocity constant");
  close(dy[4],0,1e-18L,"EdS potential constant");
  // Literal actual background acceleration forcing is retained rather than
  // invoking exact closure to delete it.
  auto defective=dust; defective.acceleration_defect=1e-7L;
  std::array<W,5> wrong{};
  irred::cosmology::detail::ideal_acoustic_derivative(y,defective,wrong);
  close(wrong[0]-dy[0],3e-7L*y[3],1e-18L,"actual background forcing");
}
void source_direction_limits() {
  // Formal radiation endpoint only: a physical photon-density change must not
  // be mistaken for an independent O(1) L/F/B defect. Lambda-map correlation
  // has its literal a^4 support, retained below in the expected derivative.
  const W a=1e-8L,G=1e-4L,x2=1e-7L,a4=a*a*a*a;
  const irred::cosmology::detail::IdealAcousticSourceCenter c{
      a,G,0,0,G,0,x2,4.L/3,4.L/3,-1,0,1.L/3};
  const auto d=irred::cosmology::detail::ideal_acoustic_source_directions(c);
  need(d.status==S::ok,"correlated source direction availability");
  auto endpoint=c;
  endpoint.a=W(.01);
  need(irred::cosmology::detail::ideal_acoustic_source_directions(endpoint).status==S::ok,
       "source directions admit the original binary64 endpoint");
  endpoint.a=W(std::nextafter(.01,std::numeric_limits<double>::infinity()));
  need(irred::cosmology::detail::ideal_acoustic_source_directions(endpoint).status!=S::ok,
       "source directions refuse the next binary64 epoch beyond support");
  close(d.gradient[0].L,-2*a4/G,1e-44L,"photon direction preserves radiation L");
  close(d.gradient[0].F,(4.L/3)*a4/G,1e-30L,"photon F derivative has only Lambda correlation");
  close(d.gradient[0].B,(4.L/3)*a4/G,1e-30L,"photon B derivative has only Lambda correlation");
  const std::array<W,5> y{x2/4,0,x2/36,-1.L/3,-2.L/3};
  const auto force=irred::cosmology::detail::ideal_acoustic_source_force(y,d.gradient[0]);
  close(force[3],2*a4/(3*G),1e-44L,"photon source velocity has no spurious radiation force");
  // Finite quotient identity, including its nonlinear remainder and an
  // invalid uncertainty interval. These controls are independent scalar
  // algebra, not an admitted whole-trajectory uncertainty result.
  const W n=3,den=7,dn=.01L,dd=-.02L,f=n/den;
  const W exact=(n+dn)/(den+dd)-f-(dn-f*dd)/den;
  const auto r=irred::cosmology::detail::ideal_acoustic_quotient_remainder(
      f,den,std::abs(dn),std::abs(dd));
  // Both the comparison quotient subtraction and its central bound are RN
  // arithmetic. Their scalar-control assembly allowance is separate from the
  // future whole-trajectory source/arithmetic shares.
  const W assembly=16*std::numeric_limits<W>::epsilon()*(1+std::abs(f));
  need(r.status==S::ok && std::abs(exact)<=r.absolute_estimate+assembly,
       "finite source quotient remainder");
  need(irred::cosmology::detail::ideal_acoustic_quotient_remainder(f,den,1,den).status!=S::ok,
       "source quotient singular interval refusal");
}
int main() {
  equation_limits();
  source_direction_limits();
  ideal_acoustic_transport_test::controls();
  auto input=fixture();
  auto prepared=prepare_ideal_acoustic(std::move(input));
  need(prepared.status()==S::ok,"actual positive photon/baryon preparation");
  need(prepared.source() && prepared.physical_mapping() && prepared.background(),"complete source owner");
  need(prepared.preparation_work().physical_mappings==1 &&
       prepared.preparation_work().background_preparations==1,"one map and background preparation");
  auto owner=std::move(prepared);
  need(prepared.status()!=S::ok && !prepared.source() && !prepared.physical_mapping() &&
       !prepared.background(),"moved-source invalidation");
  owner=std::move(owner);
  need(owner.status()==S::ok && owner.source(),"self-move preserves owner");
  failure_witness_limits(owner);
  const std::array<double,2> k{.01,1e-4};
  const std::array<double,2> a{1e-3,.01};
  IdealAcousticPolicy augmented;
  augmented.maximum_rhs_per_wavenumber=2000000;
  augmented.maximum_rhs_batch=4000000;
  const auto batch=owner.evaluate(k,a,acoustic_all_outputs,augmented);
  // Preserve the failed batch before asserting admission: its owner counters,
  // attempted source states and finite computed fields explain a refusal.
  std::cout<<std::setprecision(std::numeric_limits<W>::max_digits10)
           <<"raw_batch status="<<static_cast<int>(batch.status)
           <<" shared="<<static_cast<int>(batch.shared_dependency_status)
           <<" rhs="<<batch.evaluation_work.rhs_evaluations
           <<" P="<<batch.evaluation_work.background_evaluations
           <<" age="<<batch.evaluation_work.age_evaluations
           <<" writes="<<batch.evaluation_work.state_element_writes
           <<" diagnostics="<<batch.evaluation_work.diagnostic_evaluations<<'\n';
  if (batch.status!=S::ok) {
    for (const auto &trajectory:batch.trajectories)
      for (std::size_t r=0;r<trajectory.attempts_started;++r) {
        const auto &trial=trajectory.attempts[r];
        std::cout<<"raw_refused_attempt k_index="<<trajectory.original_k_index
                 <<" attempt="<<r<<" status="<<static_cast<int>(trial.status)
                 <<" rhs="<<trial.work.rhs_evaluations
                 <<" P="<<trial.work.background_evaluations
                 <<" writes="<<trial.work.state_element_writes
                 <<" raw_H="<<trial.maximum_direct_hamiltonian_residual
                 <<" normalized_H="<<trial.maximum_normalized_hamiltonian_residual;
        print_failure(trial); std::cout<<'\n';
      }
    for (const auto &row:batch.rows) {
      std::cout<<"raw_refused_epoch k_index="<<row.original_k_index
               <<" a_index="<<row.original_a_index
               <<" status="<<static_cast<int>(row.epoch.status)
               <<" requested_a="<<row.requested_scale_factor
               <<" actual_a="<<row.epoch.state_scale_factor
               <<" Hcal="<<row.epoch.hcal_mpc_inverse
               <<" eta="<<row.epoch.conformal_age_mpc<<'\n';
      const auto print=[](const char *name,const auto &v) {
        std::cout<<' '<<name<<'=';
        if (v) std::cout<<*v; else std::cout<<"null";
      };
      for (std::size_t f=0;f<row.outputs.size();++f) {
        const auto &v=row.outputs[f];
        std::cout<<"raw_refused_field k_index="<<row.original_k_index
                 <<" a_index="<<row.original_a_index<<" field="<<f
                 <<" status="<<static_cast<int>(v.status);
        print("computed",v.computed); print("value",v.value);
        std::cout<<" time="<<v.error.time_refinement
                 <<" initial="<<v.error.initial_refinement;
        print("source",v.error.common_source_background_age);
        print("arithmetic",v.error.arithmetic_storage_constraint);
        print("total",v.error.absolute_error_estimate); std::cout<<'\n';
      }
    }
  }
  need(batch.status==S::ok && batch.shared_dependency_status==S::ok,
       "explicit augmented source-inclusive native consumer");
  need(batch.rows.size()==4 && batch.trajectories.size()==2,"coarse grid shape");
  need(batch.evaluation_work.physical_mappings==0 && batch.evaluation_work.background_preparations==0,
       "evaluation reuses physical source");
  for (std::size_t i=0;i<k.size();++i) {
    const auto &trajectory=batch.trajectories[i];
    need(trajectory.attempts_started==5,"all original five attempts retained");
    for (std::size_t r=0;r<trajectory.attempts_started;++r) {
      const auto &trial=trajectory.attempts[r];
      std::cout<<"raw_attempt k_index="<<i<<" attempt="<<r
               <<" rhs="<<trial.work.rhs_evaluations
               <<" P="<<trial.work.background_evaluations
               <<" writes="<<trial.work.state_element_writes
               <<" raw_H="<<trial.maximum_direct_hamiltonian_residual
               <<" reduced_C="<<trial.maximum_absolute_hamiltonian_residual
               <<" normalized_H="<<trial.maximum_normalized_hamiltonian_residual
               <<" assembly="<<trial.maximum_hamiltonian_assembly_discrepancy<<'\n';
    }
    for (const auto &attempt:trajectory.attempts) {
      need(attempt.status==S::ok,"native acoustic numerical attempt");
      need(!attempt.failure,"successful attempt has no first-refusal witness");
      need(attempt.unprojected_initial_state && attempt.projected_initial_state &&
           attempt.initial_delta_projection,"initial source and actual projection retained");
    }
    for (std::size_t j=0;j<a.size();++j) {
      const auto &row=batch.rows[i*a.size()+j];
      need(row.original_k_index==i && row.original_a_index==j &&
           row.wavenumber_mpc_inverse==k[i] && row.requested_scale_factor==a[j],"caller axes preserved");
      need(row.epoch.status==S::ok && row.epoch.hcal_mpc_inverse>0 &&
           row.epoch.conformal_age_mpc>0,"matching computed epoch");
      for (const auto &value:row.outputs)
        need(value.computed && std::isfinite(*value.computed) && value.status==S::ok &&
             value.value && value.error.common_source_background_age &&
             value.error.arithmetic_storage_constraint &&
             value.error.absolute_error_estimate,"complete conditional signed prediction");
      need(row.epoch.age_error_mpc && row.epoch.derivative_consistency_estimate,
           "complete shared eta/derivative witnesses");
      std::cout<<"raw_fields k_index="<<i<<" a_index="<<j
               <<" state_a="<<row.epoch.state_scale_factor
               <<" Hcal="<<row.epoch.hcal_mpc_inverse
               <<" eta="<<row.epoch.conformal_age_mpc;
      for (const auto &value:row.outputs) std::cout<<' '<<*value.computed;
      std::cout<<'\n';
      // No missing common-source diagnostic may become an accepted zero.
      if (batch.shared_dependency_status!=S::ok)
        for (const auto &value:row.outputs)
          need(!value.value && !value.error.common_source_background_age &&
               !value.error.absolute_error_estimate,"unearned scalar admission withheld");
      close(*row.outputs[5].computed,3* *row.outputs[2].computed,2e-16L,
            "zero baryon-photon entropy identity");
    }
  }
  const std::array<double,1> one_a{.01},one_k{1e-4};
  const std::array<double,1> large_k{.01};
  const auto default_refusal=owner.evaluate(large_k,one_a,acoustic_all_outputs);
  need(default_refusal.status==S::work_limit &&
       default_refusal.trajectories[0].attempts_started>0 &&
       default_refusal.trajectories[0].attempts_started<5 &&
       default_refusal.evaluation_work.rhs_evaluations<=1000000 &&
       default_refusal.evaluation_work.rhs_evaluations+14>1000000 &&
       !default_refusal.rows[0].outputs[0].value,
       "original default RHS cap refuses augmented solve without owner refund");
  auto excessive=augmented; excessive.maximum_rhs_batch=4000001;
  need(owner.evaluate(one_k,one_a,1,excessive).status==S::invalid_input,
       "explicit augmented ceiling remains bounded");
  const auto masked=owner.evaluate(one_k,one_a,acoustic_intrinsic_temperature);
  need(masked.rows.size()==1 && masked.rows[0].outputs[2].computed,"requested intrinsic field");
  for (unsigned f=0;f<acoustic_output_count;++f)
    if (f!=2) need(!masked.rows[0].outputs[f].computed && !masked.rows[0].outputs[f].value,
                  "unrequested coordinate remains absent");
  IdealAcousticPolicy publication_limited;
  publication_limited.maximum_state_element_writes=
      masked.trajectories[0].work.state_element_writes;
  const auto publication_refusal=owner.evaluate(one_k,one_a,
      acoustic_intrinsic_temperature,publication_limited);
  need(publication_refusal.status==S::work_limit && publication_refusal.rows.size()==1 &&
       publication_refusal.rows[0].outputs[2].status==S::work_limit &&
       !publication_refusal.rows[0].outputs[2].computed &&
       !publication_refusal.rows[0].outputs[2].value,
       "publication cap refusal agrees with batch status after retained solves");
  const std::array<double,2> repeated_a{.001,.001};
  need(owner.evaluate(one_k,repeated_a,1).status==S::outside_domain,"repeated time axis refusal");
  const std::array<double,1> invalid_k{0};
  const auto refused=owner.evaluate(invalid_k,one_a,1);
  need(refused.rows.size()==1 && refused.rows[0].outputs[0].status==S::outside_domain &&
       !refused.rows[0].outputs[0].computed,"invalid k retains refused row");
  IdealAcousticPolicy limited;
  limited.maximum_state_element_writes=0;
  const auto cap=owner.evaluate(one_k,one_a,1,limited);
  need(cap.rows.size()==1 && cap.rows[0].outputs[0].status==S::work_limit &&
       !cap.rows[0].outputs[0].value,"first scalar cap refusal retained");
  auto species=fixture(); species.model.species.push_back({0,1,2});
  auto species_owner=prepare_ideal_acoustic(std::move(species));
  need(species_owner.status()==S::outside_domain && species_owner.source(),"explicit relic refusal preserves source");
  auto malformed=fixture(); malformed.source_origin.clear();
  auto malformed_owner=prepare_ideal_acoustic(std::move(malformed));
  need(malformed_owner.status()==S::invalid_input && malformed.model.tcmb_kelvin==2.7,
       "preflight preserves caller request");
  std::cout<<"ideal_acoustic_contract raw physical attempts/analytic limits/source gates\n";
}
