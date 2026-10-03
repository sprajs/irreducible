#include <irred/ideal_acoustic.hpp>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
int main() {
  using namespace irred::cosmology;
  using S=irred::numerics::Status;
  const std::array<double,2> k{1e-4,.01};
  const std::array<double,1> a{.01};
  const std::array<const char*,9> fields{"Delta_c","Delta_b","Theta0","phi",
      "delta_c","delta_b","theta_A_over_k","theta_c_over_k","phi_N"};
  IdealAcousticRequest request{{70,.02,.10,2.7,0,{}},1e-10,
      "synthetic:original-ideal-acoustic150249-grid-e8bf"};
  auto owner=prepare_ideal_acoustic(std::move(request));
  std::cout<<std::setprecision(std::numeric_limits<long double>::max_digits10);
  std::cout<<"{\"model\":\""<<ideal_acoustic_model_id<<"\",\"method\":\""
           <<ideal_acoustic_method_id<<"\",\"H0_km_s_Mpc\":70,\"omega_b_physical\":0.02,"
             "\"omega_c_physical\":0.10,\"Tcmb_K\":2.7,\"species_count\":0,"
             "\"other_massless\":0,\"ai\":1e-10,\"zeta_asymptotic\":1,\"preparation_status\":"
           <<static_cast<int>(owner.status())<<"}\n";
  if (owner.status()!=S::ok || !owner.source() || !owner.background() ||
      !owner.physical_mapping() || !owner.physical_mapping()->scalar_witnesses) return 1;
  auto retained=std::move(owner);
  if (owner.status()==S::ok || owner.source() || owner.background() || owner.physical_mapping()) return 2;
  const auto result=retained.evaluate(k,a,acoustic_all_outputs);
  const auto &work=result.evaluation_work;
  std::cout<<"{\"batch_status\":"<<static_cast<int>(result.status)
           <<",\"shared_dependency_status\":"<<static_cast<int>(result.shared_dependency_status)
           <<",\"mask\":"<<result.requested_outputs<<",\"rhs\":"<<work.rhs_evaluations
           <<",\"P_queries\":"<<work.background_evaluations<<",\"age_queries\":"<<work.age_evaluations
           <<",\"state_writes\":"<<work.state_element_writes
           <<",\"maps\":"<<result.preparation_work.physical_mappings
           <<",\"preparations\":"<<result.preparation_work.background_preparations
           <<",\"payload\":";
  if (result.peak_owned_payload_bound) std::cout<<*result.peak_owned_payload_bound;
  else std::cout<<"null";
  std::cout<<"}\n";
  bool admitted=result.status==S::ok && result.shared_dependency_status==S::ok &&
      result.rows.size()==2 && result.trajectories.size()==2;
  if (result.rows.size()!=2 || result.trajectories.size()!=2) return 3;
  for (std::size_t i=0;i<result.trajectories.size();++i) {
    const auto &trajectory=result.trajectories[i];
    if (trajectory.original_k_index!=i || trajectory.wavenumber_mpc_inverse!=k[i]) return 4;
    admitted=admitted && trajectory.attempts_started==5;
    for (std::size_t r=0;r<trajectory.attempts_started;++r) {
      const auto &trial=trajectory.attempts[r];
      std::cout<<"{\"k_index\":"<<i<<",\"attempt\":"<<r<<",\"ai\":"<<trial.initial_scale_factor
               <<",\"step\":"<<trial.maximum_log_step<<",\"status\":"<<static_cast<int>(trial.status)
               <<",\"projection\":";
      if (trial.initial_delta_projection) std::cout<<*trial.initial_delta_projection;
      else std::cout<<"null";
      const auto state=[&](const char *id,const auto &v) {
        std::cout<<",\""<<id<<"\":";
        if (!v) { std::cout<<"null"; return; }
        std::cout<<'[';
        for (unsigned j=0;j<5;++j) { if (j) std::cout<<','; std::cout<<(*v)[j]; }
        std::cout<<']';
      };
      state("unprojected",trial.unprojected_initial_state);
      state("projected",trial.projected_initial_state);
      std::cout<<"}\n";
      admitted=admitted && trial.status==S::ok && trial.unprojected_initial_state &&
               trial.projected_initial_state && trial.initial_delta_projection;
    }
  }
  for (std::size_t i=0;i<result.rows.size();++i) {
    const auto &row=result.rows[i];
    if (row.original_k_index!=i || row.original_a_index!=0 ||
        row.wavenumber_mpc_inverse!=k[i] || row.requested_scale_factor!=a[0]) return 5;
    std::cout<<"{\"k_index\":"<<i<<",\"k_Mpc_inverse\":"<<k[i]<<",\"a_requested\":"<<a[0]
             <<",\"a_state\":"<<row.epoch.state_scale_factor<<",\"Hcal_Mpc_inverse\":"
             <<row.epoch.hcal_mpc_inverse<<",\"eta_Mpc\":"<<row.epoch.conformal_age_mpc<<"}\n";
    for (unsigned f=0;f<9;++f) {
      const auto &out=row.outputs[f];
      const bool available=out.status==S::ok && out.value && out.error.absolute_error_estimate;
      std::cout<<"{\"k_index\":"<<i<<",\"field\":\""<<fields[f]<<"\",\"status\":"
               <<static_cast<int>(out.status)<<",\"computed\":";
      if (out.computed && std::isfinite(*out.computed)) std::cout<<*out.computed;
      else std::cout<<"null";
      std::cout<<",\"admitted\":"<<(available ? "true" : "false")<<",\"value\":";
      if (available) std::cout<<*out.value; else std::cout<<"null";
      std::cout<<",\"error\":";
      if (available) std::cout<<*out.error.absolute_error_estimate; else std::cout<<"null";
      std::cout<<"}\n";
      admitted=admitted && available;
    }
  }
  // This actual installed consumer refuses an unearned scientific allocation;
  // successful compilation and a finite raw witness do not make it qualified.
  return admitted ? 0 : 6;
}
