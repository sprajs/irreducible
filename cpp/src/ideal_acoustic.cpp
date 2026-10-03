#include "irred/ideal_acoustic.hpp"
#include "payload_accounting.hpp"
#include "ideal_acoustic_equations.hpp"
#include "thermal_conformal_epoch.hpp"
#include "thermal_ruler.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <utility>

namespace irred::cosmology {
namespace {
using S = numerics::Status;
using W = long double;
using State = std::array<W, 6>; // Delta,S,dV,Vc,phi,eta(Mpc).
using Epoch = detail::ThermalConformalEpoch;
bool arithmetic() noexcept {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
}
bool policy(const IdealAcousticPolicy &p) noexcept {
  return arithmetic() && std::isfinite(p.absolute_tolerance) &&
         p.absolute_tolerance > 0 && p.absolute_tolerance <= 1e-8 && std::isfinite(p.relative_tolerance) &&
         p.relative_tolerance > 0 && p.relative_tolerance <= 3e-5 && std::isfinite(p.maximum_log_step) &&
         p.maximum_log_step > 0 && p.maximum_log_step <= .02 &&
         std::isfinite(p.maximum_constraint_residual) &&
         p.maximum_constraint_residual > 0 && p.maximum_constraint_residual <= 1e-6 &&
         p.maximum_wavenumbers > 0 && p.maximum_wavenumbers <= 2 &&
         p.maximum_samples > 0 && p.maximum_samples <= 4096 &&
         p.maximum_rhs_per_wavenumber <= 1000000 &&
         p.maximum_rhs_batch <= 2000000 &&
         p.maximum_background_evaluations <= 3000000 &&
         p.maximum_age_evaluations <= 1000000 &&
         p.maximum_state_element_writes <= 30000000 &&
         p.maximum_diagnostic_evaluations <= 3000000 &&
         p.maximum_native_bytes <= 16*1024*1024;
}
IdealAcousticWork difference(const IdealAcousticWork &a,
                            const IdealAcousticWork &b) noexcept {
  return {a.physical_mappings-b.physical_mappings,
          a.background_preparations-b.background_preparations,
          a.attempts-b.attempts,a.rhs_evaluations-b.rhs_evaluations,
          a.background_evaluations-b.background_evaluations,
          a.age_evaluations-b.age_evaluations,
          a.state_element_writes-b.state_element_writes,
          a.endpoint_evaluations-b.endpoint_evaluations,
          a.diagnostic_evaluations-b.diagnostic_evaluations,
          a.output_evaluations-b.output_evaluations};
}
struct Budget {
  const IdealAcousticPolicy &p;
  IdealAcousticWork &work;
  std::size_t k_rhs_start = 0;
  bool add(std::size_t &n, std::size_t count, std::size_t cap) noexcept {
    if (n > cap || count > cap-n) return false;
    n += count; return true;
  }
  bool writes(std::size_t n) noexcept {
    return add(work.state_element_writes,n,p.maximum_state_element_writes);
  }
  bool rhs() noexcept {
    // The physical vector and eta clock are two distinct derivative owners.
    if (work.rhs_evaluations-k_rhs_start > p.maximum_rhs_per_wavenumber ||
        p.maximum_rhs_per_wavenumber-(work.rhs_evaluations-k_rhs_start) < 2)
      return false;
    return add(work.rhs_evaluations,2,p.maximum_rhs_batch);
  }
  bool background(bool age=false) noexcept {
    if (age && work.age_evaluations >= p.maximum_age_evaluations) return false;
    if (!add(work.background_evaluations,1,p.maximum_background_evaluations)) return false;
    if (age) ++work.age_evaluations;
    return true;
  }
  bool diagnostic() noexcept {
    return add(work.diagnostic_evaluations,1,p.maximum_diagnostic_evaluations);
  }
};
struct Context {
  const ThermalBackground &background;
  const detail::ThermalBaryonLoading &loading;
  Budget &budget;
  W k;
};
Epoch epoch(Context &c, W a) {
  if (!c.budget.background()) { Epoch e; e.status=S::work_limit; return e; }
  ThermalPolicy p;
  p.momentum_method=c.background.momentum_method();
  return detail::thermal_conformal_epoch(c.background,a,c.k,
          W(c.background.source().omega_gamma),p);
}
// Original MB64/67/70 ideal exchange cancellation, transformed into the
// cancellation-safe comoving variables. The actual supplied-L defect remains.
S rhs_at(Context &c, W a, const Epoch &e, const State &y, State &dy) {
  const auto sound=detail::thermal_baryon_sound(c.loading,a);
  if (sound.status!=S::ok) return sound.status;
  if (!c.budget.rhs() || !c.budget.writes(dy.size())) return S::work_limit;
  const detail::IdealAcousticCoefficients coefficients{
      e.x2,e.g,e.enthalpy_fraction,e.fb+4*e.fg/3,
      sound.loading/sound.denominator,sound.sound_speed_squared_over_c_squared,
      e.acceleration_defect};
  detail::ideal_acoustic_derivative(std::span<const W,5>(y.data(),5),coefficients,
                                   std::span<W,5>(dy.data(),5));
  dy[5]=1/e.hcal;
  for (W value:dy) if (!std::isfinite(value)) return S::overflow;
  return S::ok;
}
struct Residual {
  W absolute=0, normalized=0, relative_reduced=0, direct=0, assembly_discrepancy=0;
};
Residual residual(const State &y,const Epoch &e) {
  const W B=e.fb+4*e.fg/3;
  const std::array<W,4> terms{e.x2*y[4],1.5L*e.enthalpy_fraction*y[0],
                             -1.5L*B*y[1],4.5L*B*y[2]};
  W sum=0,reduced_scale=0;
  for (W t:terms) { sum+=t; reduced_scale+=std::abs(t); }
  const W z=-y[4]+1.5L*e.enthalpy_fraction*y[3]+1.5L*B*y[2];
  const W dc=y[0]-3*y[3],db=dc-y[1],dg=4*(y[0]-3*y[3]-y[1])/3;
  const W direct=e.x2*y[4]+3*(z+y[4])+1.5L*(e.fc*dc+e.fb*db+e.fg*dg);
  const W scale=std::abs(e.x2*y[4])+3*std::abs(z+y[4])+
                1.5L*(e.fc*std::abs(dc)+e.fb*std::abs(db)+e.fg*std::abs(dg));
  auto ratio=[&](W denominator) {
    return denominator==0 ? (sum==0 ? 0 : std::numeric_limits<W>::infinity())
                          : std::abs(sum)/denominator;
  };
  return {std::abs(sum),ratio(scale),ratio(reduced_scale),direct,std::abs(direct-sum)};
}
struct AgeContext {
  Context &c;
  W endpoint, maximum_p_relative_error=0;
  S status=S::ok;
};
double age_integrand(double u,const void *opaque) {
  auto &ac=*const_cast<AgeContext*>(static_cast<const AgeContext*>(opaque));
  if (ac.status!=S::ok) return 0;
  if (!ac.c.budget.background(true)) { ac.status=S::work_limit; return 0; }
  ThermalPolicy p; p.momentum_method=ac.c.background.momentum_method();
  const auto scaled=ac.c.background.scaled_expansion(ac.endpoint*W(u),p);
  if (scaled.status!=S::ok) { ac.status=scaled.status; return 0; }
  if (!(scaled.a4_e2>scaled.error_estimate) || !(scaled.a4_e2>0)) {
    ac.status=S::conditioning_budget_exceeded; return 0;
  }
  ac.maximum_p_relative_error=std::max(ac.maximum_p_relative_error,
                                     scaled.error_estimate/scaled.a4_e2);
  return static_cast<double>(1/std::sqrt(scaled.a4_e2));
}
struct InitialAge { S status=S::invalid_input; W value=0, estimate=0; };
InitialAge initial_age(Context &c,W start) {
  AgeContext ac{c,start};
  numerics::IntegrationPolicy p{1e-11,1e-11,
    std::min(c.budget.p.maximum_age_evaluations-c.budget.work.age_evaluations,
             c.budget.p.maximum_background_evaluations-c.budget.work.background_evaluations),30};
  const auto q=numerics::integrate(age_integrand,&ac,0,1,p);
  InitialAge out;
  out.status=ac.status==S::ok ? q.status : ac.status;
  if (out.status!=S::ok) return out;
  const W scale=detail::thermal_conformal_c_km_s/
                W(c.background.source().h0_km_s_mpc)*start;
  out.value=scale*W(q.value);
  out.estimate=scale*W(q.error_estimate)+out.value*
               (ac.maximum_p_relative_error+64*std::numeric_limits<double>::epsilon());
  if (!(out.value>0) || !(out.estimate>=0) || !std::isfinite(out.value) ||
      !std::isfinite(out.estimate)) out.status=S::outside_domain;
  return out;
}
struct Snapshot {
  S status=S::invalid_input;
  State state{};
  W a=0;
  Epoch epoch;
};
struct Run {
  S status=S::invalid_input;
  IdealAcousticAttempt attempt;
  std::vector<Snapshot> samples;
};
Run evolve(Context &c,W start,std::span<const double> targets,W step,
           const InitialAge &age) {
  Run out;
  const auto before=c.budget.work;
  ++c.budget.work.attempts;
  out.attempt.initial_scale_factor=static_cast<double>(start);
  out.attempt.maximum_log_step=static_cast<double>(step);
  auto finish=[&](S cause) {
    out.status=cause; out.attempt.status=cause;
    out.attempt.work=difference(c.budget.work,before);
  };
  if (age.status!=S::ok) { finish(age.status); return out; }
  auto current=epoch(c,start);
  if (current.status!=S::ok) { finish(current.status); return out; }
  const auto &model=c.background.source();
  const W m=(W(model.omega_b)+W(model.omega_cdm))*start/W(model.omega_gamma);
  const auto sound=detail::thermal_baryon_sound(c.loading,start);
  if (sound.status!=S::ok) { finish(sound.status); return out; }
  if (m>1e-4L || sound.loading>1e-4L || current.x2>1e-6L || current.fl>1e-10L) {
    finish(S::outside_domain); return out;
  }
  constexpr W p0=-2.L/3;
  if (!c.budget.writes(6+10)) { finish(S::work_limit); return out; }
  State y{-3*p0*current.x2/8,-p0*current.x2*current.x2/96,
          -p0*current.x2/24,
          p0*(.5L+m/16-7*m*m/160-current.x2/120),
          p0*(1-m/16+7*m*m/160-current.x2/30),age.value};
  out.attempt.unprojected_initial_state=std::array<W,5>{y[0],y[1],y[2],y[3],y[4]};
  const W old_delta=y[0],B=current.fb+4*current.fg/3;
  if (!(current.enthalpy_fraction>0)) { finish(S::outside_domain); return out; }
  if (!c.budget.writes(1+10+1)) { finish(S::work_limit); return out; }
  y[0]=(-(2.L/3)*current.x2*y[4]+B*y[1]-3*B*y[2])/current.enthalpy_fraction;
  out.attempt.projected_initial_state=std::array<W,5>{y[0],y[1],y[2],y[3],y[4]};
  out.attempt.initial_delta_projection=y[0]-old_delta;
  W n=std::log(start),a=start;
  try { out.samples.reserve(targets.size()); }
  catch (const std::bad_alloc &) { finish(S::work_limit); return out; }
  for (double target:targets) {
    const W end=std::log(W(target));
    while (n<end) {
      const W h=std::min(end-n,step/(1+std::sqrt(current.x2)));
      if (!(h>0) || n+h==n) { finish(S::conditioning_budget_exceeded); return out; }
      const W next=std::min(end,n+h),mid_a=std::exp(n+h/2),next_a=std::exp(next);
      auto middle=epoch(c,mid_a),last=epoch(c,next_a);
      if (middle.status!=S::ok || last.status!=S::ok) {
        finish(middle.status!=S::ok ? middle.status : last.status); return out;
      }
      if (!c.budget.writes(30)) { finish(S::work_limit); return out; }
      State k1{},k2{},k3{},k4{},tmp{};
      S cause=rhs_at(c,a,current,y,k1);
      if (cause!=S::ok) { finish(cause); return out; }
      if (!c.budget.writes(6)) { finish(S::work_limit); return out; }
      for (unsigned j=0;j<6;++j) tmp[j]=y[j]+h*k1[j]/2;
      cause=rhs_at(c,mid_a,middle,tmp,k2);
      if (cause!=S::ok) { finish(cause); return out; }
      if (!c.budget.writes(6)) { finish(S::work_limit); return out; }
      for (unsigned j=0;j<6;++j) tmp[j]=y[j]+h*k2[j]/2;
      cause=rhs_at(c,mid_a,middle,tmp,k3);
      if (cause!=S::ok) { finish(cause); return out; }
      if (!c.budget.writes(6)) { finish(S::work_limit); return out; }
      for (unsigned j=0;j<6;++j) tmp[j]=y[j]+h*k3[j];
      cause=rhs_at(c,next_a,last,tmp,k4);
      if (cause!=S::ok) { finish(cause); return out; }
      if (!c.budget.writes(6)) { finish(S::work_limit); return out; }
      for (unsigned j=0;j<6;++j) y[j]+=h*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;
      // Commit the epoch with the new state before any subsequent admission.
      n=next; a=next_a; current=std::move(last);
      for (W value:y) if (!std::isfinite(value)) { finish(S::overflow); return out; }
      if (!c.budget.diagnostic()) { finish(S::work_limit); return out; }
      const auto constraint=residual(y,current);
      out.attempt.maximum_normalized_hamiltonian_residual=std::max(
          out.attempt.maximum_normalized_hamiltonian_residual,constraint.normalized);
      out.attempt.maximum_absolute_hamiltonian_residual=std::max(
          out.attempt.maximum_absolute_hamiltonian_residual,constraint.absolute);
      out.attempt.maximum_relative_reduced_constraint=std::max(
          out.attempt.maximum_relative_reduced_constraint,constraint.relative_reduced);
      out.attempt.maximum_direct_hamiltonian_residual=std::max(
          out.attempt.maximum_direct_hamiltonian_residual,std::abs(constraint.direct));
      out.attempt.maximum_hamiltonian_assembly_discrepancy=std::max(
          out.attempt.maximum_hamiltonian_assembly_discrepancy,constraint.assembly_discrepancy);
      out.attempt.maximum_absolute_closure_defect=std::max(
          out.attempt.maximum_absolute_closure_defect,std::abs(current.closure_defect));
      out.attempt.maximum_absolute_acceleration_defect=std::max(
          out.attempt.maximum_absolute_acceleration_defect,std::abs(current.acceleration_defect));
      if (!std::isfinite(constraint.normalized) || constraint.normalized>c.budget.p.maximum_constraint_residual) {
        finish(S::conditioning_budget_exceeded); return out;
      }
      if (!(y[5]>=0) || c.k*(y[5]+age.estimate)>20) {
        finish(S::outside_domain); return out;
      }
    }
    ++c.budget.work.endpoint_evaluations;
    if (!c.budget.writes(12)) { finish(S::work_limit); return out; }
    out.samples.push_back({S::ok,y,a,current});
  }
  finish(S::ok); return out;
}
std::array<W,9> readout(const Snapshot &s) {
  const auto &y=s.state; const W x=std::sqrt(s.epoch.x2);
  const W theta0=(y[0]-y[1])/3-y[3];
  const W phi_n=-y[4]+1.5L*s.epoch.enthalpy_fraction*y[3]+
                1.5L*(s.epoch.fb+4*s.epoch.fg/3)*y[2];
  return {y[0],y[0]-y[1]+3*y[2],theta0,y[4],y[0]-3*y[3],
          y[0]-3*y[3]-y[1],x*(y[3]+y[2]),x*y[3],phi_n};
}
std::optional<std::size_t> retained_payload(const IdealAcousticRequest &s) noexcept {
  irred::detail::PayloadAccounting bytes(sizeof(IdealAcousticTransfer));
  bytes.vector(s.model.species); bytes.string(s.source_origin);
  return bytes.result();
}
std::optional<std::size_t> peak_payload(const IdealAcousticRequest &s,
                                      std::size_t nk,std::size_t na) noexcept {
  irred::detail::PayloadAccounting bytes(sizeof(IdealAcousticTransfer)+
                                        sizeof(IdealAcousticBatch)+16384);
  bytes.vector(s.model.species); bytes.string(s.source_origin);
  if (nk && na>std::numeric_limits<std::size_t>::max()/nk) return {};
  bytes.add(nk*na,sizeof(IdealAcousticRow));
  bytes.add(nk,sizeof(IdealAcousticTrajectory));
  // Five retained runs plus current return/snapshot lifetimes. Epoch owns its
  // raw query and optional diagnostics; charge actual sizeof(Snapshot).
  bytes.add(7*na,sizeof(Snapshot));
  bytes.add(7,sizeof(Run));
  return bytes.result();
}
} // namespace
std::optional<std::size_t> IdealAcousticTransfer::retained_payload_bound() const noexcept {
  return source_ ? retained_payload(*source_) : std::optional<std::size_t>{};
}
IdealAcousticTransfer prepare_ideal_acoustic(IdealAcousticRequest &&request,
                                             IdealAcousticPolicy p) {
  IdealAcousticTransfer out;
  const auto bytes=retained_payload(request);
  if (!policy(p) || request.source_origin.empty() || !bytes ||
      *bytes>p.maximum_native_bytes || !std::isfinite(request.initial_scale_factor) ||
      !(request.initial_scale_factor>=1e-20) || request.initial_scale_factor>1e-4)
    return out; // Structural preflight retains the original caller request.
  out.source_=std::move(request);
  const auto &m=out.source_->model;
  if (!m.species.empty() || m.physical_massless_nonphoton_density!=0 ||
      !(m.physical_baryon_density>0) || !(m.physical_cdm_density>0) ||
      !(m.tcmb_kelvin>0)) { out.status_=S::outside_domain; return out; }
  ++out.preparation_work_.physical_mappings;
  out.mapping_=map_thermal_physical_model(m);
  if (out.mapping_->status!=S::ok || !out.mapping_->model ||
      !out.mapping_->scalar_witnesses) { out.status_=out.mapping_->status; return out; }
  ++out.preparation_work_.background_preparations;
  out.background_=prepare_thermal_background(*out.mapping_->model);
  if (out.background_->status()!=S::ok) { out.status_=out.background_->status(); return out; }
  const auto coefficient=detail::ThermalRetainedCoefficientAccess::capture(*out.background_);
  const auto loading=detail::thermal_baryon_loading(*out.background_);
  if (!coefficient || loading.status!=S::ok) {
    out.status_=coefficient ? loading.status : S::conditioning_budget_exceeded; return out;
  }
  out.loading_ratio_=loading.ratio_today; out.loading_error_=loading.arithmetic_estimate;
  out.loading_numerator_=loading.numerator; out.loading_denominator_=loading.denominator;
  out.retained_lambda_=coefficient->lambda_retained;
  out.lambda_getter_signed_loss_=coefficient->lambda_getter_signed_loss;
  out.retained_critical_density_=coefficient->critical_density_ev4;
  out.retained_species_today_=coefficient->omega_species_today;
  out.retained_species_normalization_error_=coefficient->species_normalization_error;
  out.retained_lambda_emitted_=coefficient->lambda_emitted;
  out.retained_momentum_method_=coefficient->momentum_method;
  out.status_=S::ok; return out;
}
IdealAcousticTransfer::IdealAcousticTransfer(IdealAcousticTransfer &&o) noexcept {
  *this=std::move(o);
}
IdealAcousticTransfer &IdealAcousticTransfer::operator=(IdealAcousticTransfer &&o) noexcept {
  if (this!=&o) {
    status_=o.status_; source_=std::move(o.source_); mapping_=std::move(o.mapping_);
    background_=std::move(o.background_); loading_ratio_=o.loading_ratio_;
    loading_error_=o.loading_error_; loading_numerator_=o.loading_numerator_;
    loading_denominator_=o.loading_denominator_; retained_lambda_=o.retained_lambda_;
    lambda_getter_signed_loss_=o.lambda_getter_signed_loss_;
    retained_critical_density_=o.retained_critical_density_;
    retained_species_today_=o.retained_species_today_;
    retained_species_normalization_error_=o.retained_species_normalization_error_;
    retained_lambda_emitted_=o.retained_lambda_emitted_;
    retained_momentum_method_=o.retained_momentum_method_;
    preparation_work_=o.preparation_work_;
    o.source_.reset(); o.mapping_.reset(); o.background_.reset();
    o.status_=S::invalid_input; o.preparation_work_={};
  }
  return *this;
}
IdealAcousticBatch IdealAcousticTransfer::evaluate(std::span<const double> ks,
                 std::span<const double> targets,unsigned outputs,IdealAcousticPolicy p) const {
  IdealAcousticBatch out;
  out.preparation_work=preparation_work_;
  if (status_!=S::ok || !source_ || !background_) { out.status=status_; return out; }
  if (!policy(p) || !outputs || (outputs&~acoustic_all_outputs) || ks.empty() || targets.empty()) return out;
  if (ks.size()>p.maximum_wavenumbers || targets.size()>p.maximum_samples ||
      ks.size()>p.maximum_samples/targets.size()) { out.status=S::work_limit; return out; }
  const auto peak=peak_payload(*source_,ks.size(),targets.size());
  if (!peak || *peak>p.maximum_native_bytes) { out.status=S::work_limit; return out; }
  out.peak_owned_payload_bound=peak;
  for (std::size_t j=0;j<targets.size();++j) {
    if (!std::isfinite(targets[j])) { out.status=S::nonfinite_input; return out; }
    if (targets[j]<source_->initial_scale_factor || targets[j]>.01 ||
        (j && !(targets[j]>targets[j-1]))) { out.status=S::outside_domain; return out; }
  }
  out.requested_outputs=outputs;
  try {
    out.rows.resize(ks.size()*targets.size()); out.trajectories.resize(ks.size());
    Budget budget{p,out.evaluation_work};
    const detail::ThermalBaryonLoading loading{S::ok,loading_ratio_,loading_error_,
                                             loading_numerator_,loading_denominator_};
    Context initial_context{*background_,loading,budget,0};
    const W ai=source_->initial_scale_factor;
    const std::array<InitialAge,3> ages{initial_age(initial_context,ai),
                                     initial_age(initial_context,ai/2),
                                     initial_age(initial_context,ai/4)};
    bool all_ok=true;
    // The algebraic helper deliberately has no complete source/derivative
    // bundle. Preserve computed runs; refuse admission until its owner supplies
    // the narrowly qualified same-source bridge and response propagation.
    out.shared_dependency_status=S::conditioning_budget_exceeded;
    for (const auto &age:ages) if (age.status!=S::ok) {
      out.shared_dependency_status=age.status; break;
    }
    S aggregate=out.shared_dependency_status;
    for (std::size_t i=0;i<ks.size();++i) {
      auto &trajectory=out.trajectories[i];
      trajectory.original_k_index=i; trajectory.wavenumber_mpc_inverse=ks[i];
      const auto before=out.evaluation_work;
      budget.k_rhs_start=out.evaluation_work.rhs_evaluations;
      S cause=S::ok;
      if (!std::isfinite(ks[i])) cause=S::nonfinite_input;
      else if (ks[i]<1e-7 || ks[i]>.01) cause=S::outside_domain;
      Context c{*background_,loading,budget,W(ks[i])};
      std::array<Run,5> runs;
      constexpr std::array<unsigned,5> age_index{0,0,0,1,2};
      const std::array<W,5> starts{ai,ai,ai,ai/2,ai/4};
      const std::array<W,5> steps{W(p.maximum_log_step),W(p.maximum_log_step)/2,
                                  W(p.maximum_log_step)/4,W(p.maximum_log_step)/4,
                                  W(p.maximum_log_step)/4};
      if (cause==S::ok) for (unsigned r=0;r<5;++r) {
        runs[r]=evolve(c,starts[r],targets,steps[r],ages[age_index[r]]);
        trajectory.attempts[r]=runs[r].attempt; ++trajectory.attempts_started;
        if (runs[r].status!=S::ok) { cause=runs[r].status; break; }
      }
      trajectory.work=difference(out.evaluation_work,before);
      if (cause!=S::ok) aggregate=cause;
      for (std::size_t j=0;j<targets.size();++j) {
        auto &row=out.rows[i*targets.size()+j];
        row.original_k_index=i; row.original_a_index=j;
        row.wavenumber_mpc_inverse=ks[i]; row.requested_scale_factor=targets[j];
        if (cause!=S::ok) {
          for (unsigned f=0;f<9;++f) if (outputs&(1u<<f)) row.outputs[f].status=cause;
          row.epoch.status=cause; all_ok=false; continue;
        }
        std::array<std::array<W,9>,5> values;
        if (!budget.writes(45)) {
          cause=S::work_limit; all_ok=false; row.epoch.status=cause;
          for (unsigned f=0;f<9;++f) if (outputs&(1u<<f)) row.outputs[f].status=cause;
          continue;
        }
        for (unsigned r=0;r<5;++r) values[r]=readout(runs[r].samples[j]);
        out.evaluation_work.output_evaluations+=5;
        const auto &sample=runs[2].samples[j];
        row.epoch.status=S::ok; row.epoch.state_scale_factor=sample.a;
        row.epoch.hcal_mpc_inverse=sample.epoch.hcal;
        row.epoch.conformal_age_mpc=sample.state[5];
        row.epoch.endpoint_scale_factor_error=std::abs(sample.a-W(targets[j]));
        for (unsigned f=0;f<9;++f) if (outputs&(1u<<f)) {
          auto &value=row.outputs[f];
          if (!std::isfinite(values[2][f])) { value.status=S::overflow; all_ok=false; continue; }
          value.computed=values[2][f];
          value.error.time_refinement=static_cast<double>(8*std::max(
              std::abs(values[2][f]-values[1][f]),std::abs(values[1][f]-values[0][f])));
          value.error.initial_refinement=static_cast<double>(8*std::max(
              std::abs(values[2][f]-values[3][f]),std::abs(values[3][f]-values[4][f])));
          value.status=S::conditioning_budget_exceeded;
        }
        all_ok=false;
      }
    }
    out.status=all_ok ? S::ok : aggregate;
  } catch (const std::bad_alloc &) {
    out.status=S::work_limit; out.shared_dependency_status=S::work_limit;
    for (auto &row:out.rows) for (auto &value:row.outputs) {
      value.value.reset(); value.status=S::work_limit;
    }
  }
  return out;
}
} // namespace irred::cosmology
