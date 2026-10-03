#include "irred/ideal_acoustic.hpp"
#include "../src/ideal_acoustic_equations.hpp"
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
int main() {
  equation_limits();
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
  const std::array<double,2> k{.01,1e-4};
  const std::array<double,2> a{1e-3,.01};
  const auto batch=owner.evaluate(k,a,acoustic_all_outputs);
  need(batch.rows.size()==4 && batch.trajectories.size()==2,"coarse grid shape");
  need(batch.evaluation_work.physical_mappings==0 && batch.evaluation_work.background_preparations==0,
       "evaluation reuses physical source");
  std::cout<<std::setprecision(std::numeric_limits<W>::max_digits10)
           <<"raw_batch status="<<static_cast<int>(batch.status)
           <<" shared="<<static_cast<int>(batch.shared_dependency_status)
           <<" rhs="<<batch.evaluation_work.rhs_evaluations
           <<" P="<<batch.evaluation_work.background_evaluations
           <<" age="<<batch.evaluation_work.age_evaluations
           <<" writes="<<batch.evaluation_work.state_element_writes<<'\n';
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
        need(value.computed && std::isfinite(*value.computed),"actual finite signed perturbation witness");
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
  const auto masked=owner.evaluate(one_k,one_a,acoustic_intrinsic_temperature);
  need(masked.rows.size()==1 && masked.rows[0].outputs[2].computed,"requested intrinsic field");
  for (unsigned f=0;f<acoustic_output_count;++f)
    if (f!=2) need(!masked.rows[0].outputs[f].computed && !masked.rows[0].outputs[f].value,
                  "unrequested coordinate remains absent");
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
