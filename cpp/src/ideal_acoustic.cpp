#include "irred/ideal_acoustic.hpp"
#include "payload_accounting.hpp"
#include "ideal_acoustic_equations.hpp"
#include "ideal_acoustic_bridge.hpp"
#include "ideal_acoustic_initial_bounds.hpp"
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
using Response = detail::IdealAcousticResponseState;
using Frame = detail::IdealAcousticTransportFrame;
using FailureStage = IdealAcousticFailureStage;
using RadiusArithmetic = detail::ideal_acoustic_transport_internal::Arithmetic;
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
         p.maximum_rhs_per_wavenumber <= 2000000 &&
         p.maximum_rhs_batch <= 4000000 &&
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
  bool rhs(std::size_t owners=14) noexcept {
    // Nominal physical/eta plus the physical/eta owners of four signed source
    // responses, one finite-source remainder and one arithmetic radius.
    if (work.rhs_evaluations-k_rhs_start > p.maximum_rhs_per_wavenumber ||
        p.maximum_rhs_per_wavenumber-(work.rhs_evaluations-k_rhs_start) < owners)
      return false;
    return add(work.rhs_evaluations,owners,p.maximum_rhs_batch);
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
  const detail::ThermalRetainedCoefficientWitness &retained;
  const std::array<ThermalScalarMapWitness,4> &mapping;
  const detail::IdealAcousticSourceUncertainty &uncertainty;
};
Epoch epoch(Context &c, W a) {
  if (!c.budget.background()) { Epoch e; e.status=S::work_limit; return e; }
  ThermalPolicy p;
  p.momentum_method=c.background.momentum_method();
  auto e=detail::thermal_conformal_epoch(c.background,a,c.k,
          W(c.background.source().omega_gamma),p);
  if (e.status==S::ok) {
    if (!c.budget.diagnostic()) e.status=S::work_limit;
    else e.status=detail::thermal_conformal_no_species_diagnostics(
        c.background,c.retained,c.mapping,a,e);
  }
  return e;
}
detail::IdealAcousticTransportAccounting accounting(Budget &b) noexcept {
  return {&b,
      [](void *p,std::size_t n) noexcept {
        return static_cast<Budget*>(p)->writes(n);
      },
      [](void *p) noexcept { return static_cast<Budget*>(p)->diagnostic(); }};
}
detail::IdealAcousticCoefficients coefficients(const Epoch &e,
       const detail::ThermalBaryonSound &sound) noexcept {
  return {e.x2,e.g,e.enthalpy_fraction,e.fb+4*e.fg/3,
          sound.loading/sound.denominator,
          sound.sound_speed_squared_over_c_squared,e.acceleration_defect};
}
S frame(Context &c,W a,const Epoch &e,
        const detail::ThermalBaryonSound &sound,W n,W n_radius,Frame &out) {
  if (!e.shadow) return S::conditioning_budget_exceeded;
  const auto actual=coefficients(e,sound);
  const auto &m=c.background.source();
  const detail::IdealAcousticSourceCenter center{
      a,W(m.omega_gamma),W(m.omega_b),W(m.omega_cdm),
      e.shadow->p_shadow,e.shadow->p_shadow_n,actual.x2,actual.F,actual.B,
      actual.L,actual.loading_over_one_plus_loading,actual.sound_speed_squared};
  const auto directions=detail::ideal_acoustic_source_directions(center);
  const detail::IdealAcousticBridgeClock clock{S::ok,n,n_radius,true};
  return detail::ideal_acoustic_transport_frame(e,c.loading,sound,actual,center,
      directions,c.retained.lambda_retained,c.uncertainty,clock,out,
      accounting(c.budget));
}
// Original MB64/67/70 ideal exchange cancellation, transformed into the
// cancellation-safe comoving variables. The actual supplied-L defect remains.
S rhs_at(Context &c, W a, const Epoch &e, W n,W n_radius,
         const State &y,const Response &r, State &dy,Response &dr,
         FailureStage &failure_stage) {
  failure_stage=FailureStage::rk_nominal_rhs;
  const auto sound=detail::thermal_baryon_sound(c.loading,a);
  if (sound.status!=S::ok) return sound.status;
  if (!c.budget.rhs() || !c.budget.writes(dy.size())) return S::work_limit;
  const auto actual=coefficients(e,sound);
  detail::ideal_acoustic_derivative(std::span<const W,5>(y.data(),5),actual,
                                   std::span<W,5>(dy.data(),5));
  dy[5]=1/e.hcal;
  for (W value:dy) if (!std::isfinite(value)) return S::overflow;
  Frame f;
  failure_stage=FailureStage::rk_frame;
  auto status=frame(c,a,e,sound,n,n_radius,f);
  if (status!=S::ok) return status;
  failure_stage=FailureStage::rk_response_rhs;
  return detail::ideal_acoustic_transport_rhs(y,r,f,c.uncertainty,dr,
                                             accounting(c.budget));
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
struct InitialAge {
  S status=S::invalid_input;
  W value=0, estimate=0, source_estimate=0;
};
InitialAge initial_age(Context &c,W start) {
  const auto &limit=c.budget.p;
  const auto &spent=c.budget.work;
  if (spent.age_evaluations>limit.maximum_age_evaluations ||
      spent.background_evaluations>limit.maximum_background_evaluations)
    return {S::work_limit,0,0};
  const auto available=std::min(
      limit.maximum_age_evaluations-spent.age_evaluations,
      limit.maximum_background_evaluations-spent.background_evaluations);
  // The quadrature requires its three initial callbacks. An exhausted caller
  // budget is an owned resource refusal, before the generic integrator sees an
  // invalid zero-callback policy; no callback or counter is silently reset.
  if (available<3) return {S::work_limit,0,0};
  AgeContext ac{c,start};
  numerics::IntegrationPolicy p{1e-11,1e-11,available,30};
  const auto q=numerics::integrate(age_integrand,&ac,0,1,p);
  InitialAge out;
  out.status=ac.status==S::ok ? q.status : ac.status;
  if (out.status!=S::ok) return out;
  const W scale=detail::thermal_conformal_c_km_s/
                W(c.background.source().h0_km_s_mpc)*start;
  out.value=scale*W(q.value);
  out.estimate=scale*W(q.error_estimate)+out.value*
               (ac.maximum_p_relative_error+64*std::numeric_limits<double>::epsilon());
  RadiusArithmetic radius;
  std::array<W,4> amplitude;
  for (unsigned j=0;j<4;++j)
    amplitude[j]=radius.plus(radius.magnitude(c.uncertainty.signed_shift[j]),
                             c.uncertainty.operation_radius[j]);
  const W a2=radius.mul(start,start),a4=radius.mul(a2,a2);
  const W vacuum=radius.plus(radius.plus(amplitude[0],amplitude[1]),
                             radius.plus(amplitude[2],amplitude[3]));
  const W dp=radius.plus(amplitude[0],radius.plus(
      radius.times(start,radius.plus(amplitude[1],amplitude[2])),
      radius.times(a4,vacuum)));
  // For every u in[0,1], nominal P>=Gamma and |deltaP|<=dp. The full
  // positive-family inverse-sqrt change is <=t/(1-t), without another H.
  const W photon=W(c.background.source().omega_gamma);
  const W lower=radius.lower(photon,dp);
  out.source_estimate=radius.times(radius.plus(out.value,out.estimate),
                                    radius.over(dp,lower));
  if (radius.status!=S::ok) { out.status=radius.status; return out; }
  if (!(out.value>0) || !(out.estimate>=0) || !std::isfinite(out.value) ||
      !std::isfinite(out.estimate)) out.status=S::outside_domain;
  return out;
}
struct Snapshot {
  S status=S::invalid_input;
  State state{};
  Response response{};
  W n=0,n_radius=0;
  W a=0;
  Epoch epoch;
};
struct Run {
  S status=S::invalid_input;
  IdealAcousticAttempt attempt;
  std::vector<Snapshot> samples;
};
// Conditional normal-Wide RK assembly. The stored positive comparison vectors
// are rounded upward; the nominal and four signed response assembly losses
// enter ARITHMETIC once. The mesh estimator remains a separate empirical gate.
S combine(Context &c,const State &base,const Response &base_r,
          const std::array<const State*,4> &slopes,
          const std::array<const Response*,4> &response_slopes,
          const std::array<W,4> &weights,State &out,Response &out_r,
          std::optional<IdealAcousticRadiusFailure> &failure) {
  if (!c.budget.writes(42) || !c.budget.diagnostic()) return S::work_limit;
  RadiusArithmetic a;
  for (unsigned i=0;i<6;++i) {
    W value=base[i],scale=std::abs(base[i]);
    for (unsigned s=0;s<4;++s) if (weights[s]!=0) {
      const W term=a.mul(weights[s],(*slopes[s])[i]);
      value=a.add(value,term); scale=a.plus(scale,a.magnitude(term));
    }
    out[i]=value;
    W assembly=a.loss(scale);
    for (unsigned j=0;j<4;++j) {
      const unsigned n=6*j+i;
      W z=base_r[n],zs=a.magnitude(z);
      for (unsigned s=0;s<4;++s) if (weights[s]!=0) {
        const W term=a.mul(weights[s],(*response_slopes[s])[n]);
        z=a.add(z,term); zs=a.plus(zs,a.magnitude(term));
      }
      out_r[n]=z;
      const W amplitude=a.plus(a.magnitude(c.uncertainty.signed_shift[j]),
                                c.uncertainty.operation_radius[j]);
      assembly=a.plus(assembly,a.times(amplitude,a.loss(zs)));
    }
    for (unsigned offset:{24u,30u}) {
      W radius=base_r[offset+i],rs=a.magnitude(radius);
      for (unsigned s=0;s<4;++s) if (weights[s]!=0) {
        const W term=a.mul(weights[s],(*response_slopes[s])[offset+i]);
        radius=a.add(radius,term); rs=a.plus(rs,a.magnitude(term));
      }
      const W provisional=a.signed_upper(a.add(radius,a.loss(rs)));
      radius=offset==30 ? a.assemble_signed_radius(provisional,assembly) : provisional;
      if (!(radius>=0) || !detail::ideal_acoustic_transport_internal::radius(radius)) {
        failure=detail::ideal_acoustic_transport_internal::record_radius_failure(
            i,static_cast<IdealAcousticRadiusChannel>(offset),provisional,
            offset==30 ? std::optional<W>{assembly} : std::nullopt,radius,a.status);
        return S::conditioning_budget_exceeded;
      }
      out_r[offset+i]=radius;
    }
  }
  return a.status;
}
detail::IdealAcousticSourceUncertainty source_uncertainty(
    const ThermalBackground &background,
    const std::array<ThermalScalarMapWitness,4> &map,
    const detail::ThermalRetainedCoefficientWitness &retained) {
  detail::IdealAcousticSourceUncertainty out;
  const auto &model=background.source();
  const std::array<double,4> emitted{model.omega_gamma,model.omega_b,
                                    model.omega_cdm,model.omega_massless_nonphoton};
  for (unsigned j=0;j<4;++j) {
    const auto &w=map[j];
    if (!std::isfinite(w.wide_value) || !std::isfinite(w.wide_operation_estimate) ||
        !std::isfinite(w.measured_absolute_cast_loss) ||
        w.emitted_value!=emitted[j] ||
        w.measured_absolute_cast_loss!=std::abs(w.wide_value-W(w.emitted_value)) ||
        w.wide_operation_estimate<0 ||
        (j<3 && (!(w.wide_value>0) || !(w.wide_operation_estimate>0))))
      return out;
    if (j<3) {
      out.signed_shift[j]=w.wide_value-W(w.emitted_value);
      out.operation_radius[j]=w.wide_operation_estimate;
    }
  }
  out.original_other_shift=map[3].wide_value-W(map[3].emitted_value);
  out.original_other_radius=map[3].wide_operation_estimate;
  out.original_other_zero_witness=map[3].wide_value==0 &&
      map[3].emitted_value==0 && map[3].wide_operation_estimate==0 &&
      map[3].measured_absolute_cast_loss==0;
  if (!out.original_other_zero_witness || retained.status!=S::ok ||
      retained.omega_species_today!=0 || retained.species_normalization_error!=0)
    return out;
  const W photon=W(model.omega_gamma),baryon=W(model.omega_b),cdm=W(model.omega_cdm);
  out.signed_shift[3]=(1-photon-baryon-cdm)-retained.lambda_retained;
  out.operation_radius[3]=32*std::numeric_limits<W>::epsilon()*
      (1+photon+baryon+cdm+std::abs(retained.lambda_retained));
  if (!std::isfinite(out.signed_shift[3]) ||
      !(out.operation_radius[3]>0) || !std::isfinite(out.operation_radius[3]))
    return out;
  out.status=S::ok; return out;
}
Run evolve(Context &c,W start,std::span<const double> targets,W step,
           const InitialAge &age) {
  Run out;
  const auto before=c.budget.work;
  ++c.budget.work.attempts;
  out.attempt.initial_scale_factor=static_cast<double>(start);
  out.attempt.maximum_log_step=static_cast<double>(step);
  std::optional<W> committed_a;
  std::optional<IdealAcousticRadiusFailure> radius_failure;
  auto finish=[&](S cause,FailureStage stage,unsigned rk_stage=0,
                  std::optional<W> coefficient_a={}) {
    out.status=cause; out.attempt.status=cause;
    out.attempt.work=difference(c.budget.work,before);
    if (cause!=S::ok)
      out.attempt.failure=IdealAcousticFailure{
          stage,rk_stage,committed_a,coefficient_a,std::move(radius_failure)};
  };
  if (age.status!=S::ok) { finish(age.status,FailureStage::initial_age); return out; }
  auto current=epoch(c,start);
  if (current.status!=S::ok) { finish(current.status,FailureStage::initial_background,0,start); return out; }
  const auto &model=c.background.source();
  const W m=(W(model.omega_b)+W(model.omega_cdm))*start/W(model.omega_gamma);
  const auto sound=detail::thermal_baryon_sound(c.loading,start);
  if (sound.status!=S::ok) { finish(sound.status,FailureStage::initial_loading,0,start); return out; }
  if (m>1e-4L || sound.loading>1e-4L || current.x2>1e-6L || current.fl>1e-10L) {
    finish(S::outside_domain,FailureStage::initial_domain,0,start); return out;
  }
  constexpr W p0=-2.L/3;
  if (!c.budget.writes(6+10)) { finish(S::work_limit,FailureStage::initial_state_storage,0,start); return out; }
  State y{-3*p0*current.x2/8,-p0*current.x2*current.x2/96,
          -p0*current.x2/24,
          p0*(.5L+m/16-7*m*m/160-current.x2/120),
          p0*(1-m/16+7*m*m/160-current.x2/30),age.value};
  out.attempt.unprojected_initial_state=std::array<W,5>{y[0],y[1],y[2],y[3],y[4]};
  const W old_delta=y[0],B=current.fb+4*current.fg/3;
  if (!(current.enthalpy_fraction>0)) { finish(S::outside_domain,FailureStage::initial_projection,0,start); return out; }
  if (!c.budget.writes(1+10+1)) { finish(S::work_limit,FailureStage::initial_projection,0,start); return out; }
  y[0]=(-(2.L/3)*current.x2*y[4]+B*y[1]-3*B*y[2])/current.enthalpy_fraction;
  out.attempt.projected_initial_state=std::array<W,5>{y[0],y[1],y[2],y[3],y[4]};
  out.attempt.initial_delta_projection=y[0]-old_delta;
  committed_a=start;
  W n=std::log(start),a=start;
  RadiusArithmetic clock_arithmetic;
  // The selected log law is relative to the true log. Outward128eps covers
  // its implicit stored-versus-true denominator under64eps<1/2.
  W n_radius=clock_arithmetic.times(2,
      clock_arithmetic.loss(clock_arithmetic.plus(1,std::abs(n))));
  if (!c.budget.writes(36)) { finish(S::work_limit,FailureStage::initial_clock,0,start); return out; }
  Response response{};
  Frame initial_frame;
  auto initial_status=frame(c,a,current,sound,n,n_radius,initial_frame);
  FailureStage initial_stage=FailureStage::initial_frame;
  detail::IdealAcousticTransportInitialBounds initial_bounds;
  if (initial_status==S::ok) {
    initial_stage=FailureStage::initial_source_bounds;
    initial_status=detail::ideal_acoustic_initial_bounds(y,initial_frame,
        c.uncertainty,c.retained.lambda_retained,age.source_estimate,age.estimate,
        initial_bounds,accounting(c.budget));
  }
  if (initial_status==S::ok) {
    initial_stage=FailureStage::initial_response;
    initial_status=detail::ideal_acoustic_transport_initial(y,initial_frame,
        c.uncertainty,initial_bounds,response,accounting(c.budget));
  }
  if (initial_status!=S::ok) { finish(initial_status,initial_stage,0,start); return out; }
  try { out.samples.reserve(targets.size()); }
  catch (const std::bad_alloc &) { finish(S::work_limit,FailureStage::sample_storage); return out; }
  for (double target:targets) {
    const W end=std::log(W(target));
    while (n<end) {
      const W h=std::min(end-n,step/(1+std::sqrt(current.x2)));
      if (!(h>0) || n+h==n) { finish(S::conditioning_budget_exceeded,FailureStage::step_clock); return out; }
      const W mid_n=n+h/2,next=std::min(end,n+h);
      const W mid_a=std::exp(mid_n),next_a=std::exp(next);
      const W mid_radius=clock_arithmetic.plus(n_radius,
          clock_arithmetic.loss(clock_arithmetic.plus(std::abs(n),h/2)));
      const W next_radius=std::max(clock_arithmetic.plus(n_radius,
          clock_arithmetic.loss(clock_arithmetic.plus(std::abs(n),h))),
          clock_arithmetic.loss(clock_arithmetic.plus(1,std::abs(next))));
      if (clock_arithmetic.status!=S::ok) {
        finish(clock_arithmetic.status,FailureStage::step_clock); return out;
      }
      auto middle=epoch(c,mid_a),last=epoch(c,next_a);
      if (middle.status!=S::ok || last.status!=S::ok) {
        finish(middle.status!=S::ok ? middle.status : last.status,
               FailureStage::step_background,middle.status!=S::ok ? 2 : 4,
               middle.status!=S::ok ? mid_a : next_a); return out;
      }
      if (!c.budget.writes(210)) { finish(S::work_limit,FailureStage::rk_storage); return out; }
      State k1{},k2{},k3{},k4{},tmp{};
      Response r1{},r2{},r3{},r4{},tmp_r{};
      FailureStage rhs_stage=FailureStage::rk_nominal_rhs;
      S cause=rhs_at(c,a,current,n,n_radius,y,response,k1,r1,rhs_stage);
      if (cause!=S::ok) { finish(cause,rhs_stage,1,a); return out; }
      cause=combine(c,y,response,{&k1,&k1,&k1,&k1},{&r1,&r1,&r1,&r1},
                    {h/2,0,0,0},tmp,tmp_r,radius_failure);
      if (cause!=S::ok) { finish(cause,FailureStage::rk_combine,2,mid_a); return out; }
      cause=rhs_at(c,mid_a,middle,mid_n,mid_radius,tmp,tmp_r,k2,r2,rhs_stage);
      if (cause!=S::ok) { finish(cause,rhs_stage,2,mid_a); return out; }
      cause=combine(c,y,response,{&k2,&k2,&k2,&k2},{&r2,&r2,&r2,&r2},
                    {h/2,0,0,0},tmp,tmp_r,radius_failure);
      if (cause!=S::ok) { finish(cause,FailureStage::rk_combine,3,mid_a); return out; }
      cause=rhs_at(c,mid_a,middle,mid_n,mid_radius,tmp,tmp_r,k3,r3,rhs_stage);
      if (cause!=S::ok) { finish(cause,rhs_stage,3,mid_a); return out; }
      cause=combine(c,y,response,{&k3,&k3,&k3,&k3},{&r3,&r3,&r3,&r3},
                    {h,0,0,0},tmp,tmp_r,radius_failure);
      if (cause!=S::ok) { finish(cause,FailureStage::rk_combine,4,next_a); return out; }
      cause=rhs_at(c,next_a,last,next,next_radius,tmp,tmp_r,k4,r4,rhs_stage);
      if (cause!=S::ok) { finish(cause,rhs_stage,4,next_a); return out; }
      if (!c.budget.writes(42)) { finish(S::work_limit,FailureStage::endpoint_combine,0,next_a); return out; }
      const auto old_y=y; const auto old_response=response;
      cause=combine(c,old_y,old_response,{&k1,&k2,&k3,&k4},{&r1,&r2,&r3,&r4},
                    {h/6,h/3,h/3,h/6},y,response,radius_failure);
      if (cause!=S::ok) { finish(cause,FailureStage::endpoint_combine,0,next_a); return out; }
      // Commit the epoch with the new state before any subsequent admission.
      n=next; n_radius=next_radius; a=next_a; current=std::move(last);
      committed_a=a;
      for (W value:y) if (!std::isfinite(value)) { finish(S::overflow,FailureStage::committed_state,0,a); return out; }
      if (!c.budget.diagnostic()) { finish(S::work_limit,FailureStage::hamiltonian_constraint,0,a); return out; }
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
        finish(S::conditioning_budget_exceeded,FailureStage::hamiltonian_constraint,0,a); return out;
      }
      if (!c.budget.diagnostic()) { finish(S::work_limit,FailureStage::phase_response,0,a); return out; }
      RadiusArithmetic phase_arithmetic;
      std::array<W,6> phase_linear{},phase_total{};
      const auto phase_status=detail::ideal_acoustic_transport_internal::response_sizes(
          response,c.uncertainty,phase_linear,phase_total,phase_arithmetic);
      if (phase_status!=S::ok) { finish(phase_status,FailureStage::phase_response,0,a); return out; }
      // Stage source/arithmetic support includes the initial age. The complete
      // mesh/start phase allowance is additionally enforced at every output.
      if (!(y[5]>=0) || c.k*phase_arithmetic.plus(y[5],phase_total[5])>20) {
        finish(S::outside_domain,FailureStage::phase_domain,0,a); return out;
      }
    }
    ++c.budget.work.endpoint_evaluations;
    if (!c.budget.writes(84)) { finish(S::work_limit,FailureStage::endpoint_storage,0,a); return out; }
    out.samples.push_back({S::ok,y,response,n,n_radius,a,current});
  }
  finish(S::ok,FailureStage::endpoint_storage); return out;
}
struct Readout {
  S status=S::invalid_input;
  std::array<W,9> value{};
  std::array<W,10> source{},arithmetic{};
  W derivative_consistency=0;
};
Readout readout(Context &c,const Snapshot &s,double target) {
  Readout out;
  const auto &y=s.state; const W x=std::sqrt(s.epoch.x2);
  const W theta0=(y[0]-y[1])/3-y[3];
  const W phi_n=-y[4]+1.5L*s.epoch.enthalpy_fraction*y[3]+
                1.5L*(s.epoch.fb+4*s.epoch.fg/3)*y[2];
  out.value={y[0],y[0]-y[1]+3*y[2],theta0,y[4],y[0]-3*y[3],
             y[0]-3*y[3]-y[1],x*(y[3]+y[2]),x*y[3],phi_n};
  // The accumulated path may differ from true log(target) even when the
  // stored exp/log roundtrip equals target. Every run lands on its stored
  // target log; own target-log error plus all accumulated path ambiguity.
  RadiusArithmetic endpoint_arithmetic;
  const W target_log_radius=endpoint_arithmetic.times(2,
      endpoint_arithmetic.loss(endpoint_arithmetic.plus(1,std::abs(s.n))));
  const W log_distance=endpoint_arithmetic.plus(s.n_radius,target_log_radius);
  if (endpoint_arithmetic.status!=S::ok) {
    out.status=endpoint_arithmetic.status; return out;
  }
  (void)target; // Original requested-a and actual-a discrepancy stay in row.
  const auto sound=detail::thermal_baryon_sound(c.loading,s.a);
  if (sound.status!=S::ok) { out.status=sound.status; return out; }
  Frame f;
  out.status=frame(c,s.a,s.epoch,sound,s.n,log_distance,f);
  if (out.status!=S::ok) return out;
  f.x=x;
  if (!c.budget.writes(36)) { out.status=S::work_limit; return out; }
  Response fixed_a=s.response;
  out.status=detail::ideal_acoustic_transport_fixed_a(y,fixed_a,f,c.uncertainty,
      log_distance,accounting(c.budget));
  if (out.status!=S::ok) return out;
  out.status=detail::ideal_acoustic_transport_readout(y,fixed_a,f,c.uncertainty,
      out.source,out.arithmetic,accounting(c.budget));
  if (out.status==S::ok && s.epoch.shadow)
    out.derivative_consistency=f.arithmetic.coefficient_radius[3];
  return out;
}
W refinement(const std::array<W,5> &v,bool initial) noexcept {
  return 8*std::max(initial ? std::abs(v[2]-v[3]) : std::abs(v[2]-v[1]),
                    initial ? std::abs(v[3]-v[4]) : std::abs(v[1]-v[0]));
}
std::optional<double> outward(W x) noexcept {
  if (!(x>=0) || !std::isfinite(x)) return {};
  double d=static_cast<double>(x);
  if (!std::isfinite(d)) return {};
  if (W(d)<x) d=std::nextafter(d,std::numeric_limits<double>::infinity());
  return std::isfinite(d) ? std::optional<double>(d) : std::nullopt;
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
    const detail::ThermalRetainedCoefficientWitness retained{
        S::ok,retained_critical_density_,retained_species_today_,retained_lambda_,
        retained_species_normalization_error_,retained_lambda_emitted_,
        lambda_getter_signed_loss_,retained_momentum_method_};
    if (!mapping_ || !mapping_->scalar_witnesses) {
      out.status=S::conditioning_budget_exceeded; return out;
    }
    const auto &map=*mapping_->scalar_witnesses;
    const auto uncertainty=source_uncertainty(*background_,map,retained);
    if (uncertainty.status!=S::ok) {
      out.status=uncertainty.status; out.shared_dependency_status=uncertainty.status;
      return out;
    }
    Context initial_context{*background_,loading,budget,0,retained,map,uncertainty};
    const W ai=source_->initial_scale_factor;
    const std::array<InitialAge,3> ages{initial_age(initial_context,ai),
                                     initial_age(initial_context,ai/2),
                                     initial_age(initial_context,ai/4)};
    bool all_ok=true;
    out.shared_dependency_status=S::ok;
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
      Context c{*background_,loading,budget,W(ks[i]),retained,map,uncertainty};
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
        std::array<Readout,5> values;
        if (!budget.writes(45)) {
          cause=S::work_limit; aggregate=cause; all_ok=false; row.epoch.status=cause;
          for (unsigned f=0;f<9;++f) if (outputs&(1u<<f)) row.outputs[f].status=cause;
          continue;
        }
        for (unsigned r=0;r<5;++r) values[r]=readout(c,runs[r].samples[j],targets[j]);
        out.evaluation_work.output_evaluations+=5;
        const auto &sample=runs[2].samples[j];
        row.epoch.status=S::ok; row.epoch.state_scale_factor=sample.a;
        row.epoch.hcal_mpc_inverse=sample.epoch.hcal;
        row.epoch.conformal_age_mpc=sample.state[5];
        row.epoch.endpoint_scale_factor_error=std::abs(sample.a-W(targets[j]));
        S response_status=S::ok;
        for (const auto &v:values) if (v.status!=S::ok) {
          response_status=v.status; break;
        }
        if (response_status==S::ok) {
          std::array<W,5> eta{},source_age{},arithmetic_age{};
          for (unsigned r=0;r<5;++r) {
            eta[r]=runs[r].samples[j].state[5];
            source_age[r]=values[r].source[9];
            arithmetic_age[r]=values[r].arithmetic[9];
          }
          RadiusArithmetic radius;
          const W source=radius.plus(source_age[2],radius.plus(
              refinement(source_age,false),refinement(source_age,true)));
          const W arithmetic=radius.plus(arithmetic_age[2],radius.plus(
              refinement(arithmetic_age,false),refinement(arithmetic_age,true)));
          const W complete_age=radius.plus(radius.plus(source,arithmetic),
              radius.plus(refinement(eta,false),refinement(eta,true)));
          if (radius.status!=S::ok) response_status=radius.status;
          else if (!(sample.state[5]>=0) || c.k*radius.plus(sample.state[5],complete_age)>20)
            response_status=S::outside_domain;
          else {
            row.epoch.age_error_mpc=complete_age;
            row.epoch.derivative_consistency_estimate=values[2].derivative_consistency;
          }
        }
        if (response_status!=S::ok) {
          row.epoch.status=response_status; aggregate=response_status; all_ok=false;
        }
        for (unsigned f=0;f<9;++f) if (outputs&(1u<<f)) {
          auto &value=row.outputs[f];
          if (!std::isfinite(values[2].value[f])) {
            value.status=S::overflow; aggregate=value.status; all_ok=false; continue;
          }
          value.computed=values[2].value[f];
          std::array<W,5> central{},source{},arithmetic{};
          for (unsigned r=0;r<5;++r) {
            central[r]=values[r].value[f]; source[r]=values[r].source[f];
            arithmetic[r]=values[r].arithmetic[f];
          }
          const auto time=outward(refinement(central,false));
          const auto initial=outward(refinement(central,true));
          if (time) value.error.time_refinement=*time;
          if (initial) value.error.initial_refinement=*initial;
          if (response_status!=S::ok || !time || !initial) {
            value.status=response_status==S::ok ? S::overflow : response_status;
            aggregate=value.status; all_ok=false; continue;
          }
          RadiusArithmetic radius;
          const W source_error=radius.plus(source[2],radius.plus(
              refinement(source,false),refinement(source,true)));
          const double stored=static_cast<double>(central[2]);
          const W cast=std::abs(central[2]-W(stored));
          const W arithmetic_error=radius.plus(cast,radius.plus(arithmetic[2],
              radius.plus(refinement(arithmetic,false),refinement(arithmetic,true))));
          const auto source_stored=outward(source_error);
          const auto arithmetic_stored=outward(arithmetic_error);
          if (!source_stored || !arithmetic_stored || !std::isfinite(stored) ||
              radius.status!=S::ok) {
            value.status=S::overflow; aggregate=value.status; all_ok=false; continue;
          }
          value.error.common_source_background_age=*source_stored;
          value.error.arithmetic_storage_constraint=*arithmetic_stored;
          // Conservative total is assembled from the STORED outward leaves.
          const W total=radius.plus(radius.plus(W(*time),W(*initial)),
              radius.plus(W(*source_stored),W(*arithmetic_stored)));
          value.error.absolute_error_estimate=outward(total);
          const W epsilon=W(p.absolute_tolerance)+W(p.relative_tolerance)*std::abs(central[2]);
          if (!value.error.absolute_error_estimate || *time>epsilon/5 ||
              *initial>epsilon/5 || *source_stored>epsilon/5 ||
              *arithmetic_stored>epsilon/5 ||
              W(*value.error.absolute_error_estimate)>4*epsilon/5) {
            value.status=S::conditioning_budget_exceeded;
            aggregate=value.status; all_ok=false; continue;
          }
          value.value=stored; value.status=S::ok;
        }
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
