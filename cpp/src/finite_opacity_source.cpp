#include "irred/finite_opacity_source.hpp"
#include "finite_opacity_transport.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace irred::cosmology {
namespace {
namespace f = detail::finite_opacity;
using W=long double; using S=numerics::Status;
bool profile() {
  return std::numeric_limits<W>::digits>=64 &&
      std::numeric_limits<W>::max_exponent>=16384 && std::fegetround()==FE_TONEAREST;
}
bool policy_valid(const FiniteOpacityPolicy &p) {
  return std::isfinite(p.state_absolute_tolerance) && p.state_absolute_tolerance>0 &&
      std::isfinite(p.relative_tolerance) && p.relative_tolerance>0 &&
      std::isfinite(p.channel_absolute_tolerance_mpc_inverse) && p.channel_absolute_tolerance_mpc_inverse>0 &&
      std::isfinite(p.maximum_eta_step_mpc) && p.maximum_eta_step_mpc>0 && p.maximum_eta_step_mpc<=4 &&
      std::isfinite(p.maximum_log_a_step) && p.maximum_log_a_step>0 && p.maximum_log_a_step<=.03 &&
      p.maximum_wavenumbers>0 && p.maximum_wavenumbers<=2 &&
      p.maximum_original_times>=2 && p.maximum_original_times<=1024 &&
      p.maximum_fine_times>=5 && p.maximum_fine_times<=4096 &&
      p.maximum_hierarchy==192 && p.maximum_attempted_steps>0 && p.maximum_attempted_steps<=100000 &&
      p.maximum_background_clock_calls>0 && p.maximum_background_clock_calls<=1000000 &&
      p.maximum_coupled_stage_solves>0 && p.maximum_coupled_stage_solves<=200000 &&
      p.maximum_destination_writes>0 && p.maximum_destination_writes<=2000000000 &&
      p.maximum_native_bytes>0 && p.maximum_native_bytes<=64*1024*1024;
}
void request_payload(irred::detail::PayloadAccounting &a,const FiniteOpacityRequest &r) {
  a.vector(r.model.species); a.vector(r.opacity.eta_mpc);
  a.vector(r.opacity.differential_opacity_mpc_inverse); a.vector(r.k_mpc_inverse);
  a.string(r.opacity.origin); a.string(r.opacity.law_id);
  a.string(r.opacity.exact_member_digest); a.string(r.source_origin);
}
S request_valid(const FiniteOpacityRequest &r,const FiniteOpacityPolicy &p) {
  if(!profile()) return S::outside_domain;
  if(!policy_valid(p)) return S::invalid_input;
  irred::detail::PayloadAccounting bytes(sizeof(FiniteOpacityIdentity)); request_payload(bytes,r);
  if(!bytes.result()) return S::overflow;
  if(*bytes.result()>p.maximum_native_bytes) return S::work_limit;
  const auto &o=r.opacity; const auto &m=r.model;
  if(!m.species.empty() || m.physical_massless_nonphoton_density!=0) return S::outside_domain;
  if(!(r.initial_scale_factor>=0x1p-512 && r.initial_scale_factor<r.observer_scale_factor &&
       r.observer_scale_factor<=1) || !std::isfinite(r.initial_scale_factor) ||
       !std::isfinite(r.observer_scale_factor)) return S::outside_domain;
  if(o.eta_mpc.size()<2 || o.eta_mpc.size()!=o.differential_opacity_mpc_inverse.size() ||
     r.k_mpc_inverse.empty() || o.origin.empty() || o.law_id.empty() ||
     o.exact_member_digest.empty() || r.source_origin.empty()) return S::invalid_input;
  if(o.eta_mpc.size()>p.maximum_original_times || r.k_mpc_inverse.size()>p.maximum_wavenumbers ||
     4*(o.eta_mpc.size()-1)+1>p.maximum_fine_times) return S::work_limit;
  // Preparation simultaneously owns its identity, scalar map/capture/epochs,
  // and requested seed/boundary storage. Account those before acquisition.
  bytes.add(sizeof(FiniteOpacitySourceProducer)+sizeof(ThermalPhysicalMapping)+
            sizeof(detail::ThermalRetainedCoefficientWitness)+2*sizeof(detail::ThermalConformalEpoch),1);
  bytes.add(r.k_mpc_inverse.size(),sizeof(FiniteOpacitySeed)+3*sizeof(W));
  if(!bytes.result())return S::overflow;
  if(*bytes.result()>p.maximum_native_bytes)return S::work_limit;
  if(!(m.physical_baryon_density>0 && m.physical_cdm_density>0 && m.tcmb_kelvin>0 && m.h0_km_s_mpc>0))
    return S::outside_domain;
  W tau=0;
  for(std::size_t i=0;i<o.eta_mpc.size();++i) {
    if(!std::isfinite(o.eta_mpc[i]) || !std::isfinite(o.differential_opacity_mpc_inverse[i])) return S::nonfinite_input;
    if(!(o.eta_mpc[i]>0) || (i && !(o.eta_mpc[i]>o.eta_mpc[i-1])) ||
       !(o.differential_opacity_mpc_inverse[i]>=0 && o.differential_opacity_mpc_inverse[i]<=1e12)) return S::outside_domain;
    if(i) tau+=(W(o.eta_mpc[i])-o.eta_mpc[i-1])*(W(o.differential_opacity_mpc_inverse[i])+o.differential_opacity_mpc_inverse[i-1])/2;
  }
  if(!std::isfinite(tau) || tau>512) return S::outside_domain;
  for(std::size_t i=0;i<r.k_mpc_inverse.size();++i) {
    const double k=r.k_mpc_inverse[i];
    if(!std::isfinite(k)) return S::nonfinite_input;
    if(!(k>0) || (i && !(k>r.k_mpc_inverse[i-1])) ||
       W(k)*(W(o.eta_mpc.back())-o.eta_mpc.front())>64) return S::outside_domain;
  }
  return S::ok;
}
detail::ThermalRetainedCoefficientWitness retained(const FiniteOpacityIdentity &id) {
  detail::ThermalRetainedCoefficientWitness r;
  r.status=S::ok; r.critical_density_ev4=id.critical_density_ev4;
  r.lambda_retained=id.lambda_retained; r.lambda_emitted=id.lambda_emitted;
  r.lambda_getter_signed_loss=id.lambda_getter_signed_loss;
  r.species_normalization_error=id.species_normalization_error;
  r.momentum_method=id.background.momentum_method(); return r;
}
detail::ThermalConformalEpoch epoch(const FiniteOpacityIdentity &id,W a,W k,f::Ledger &ledger) {
  detail::ThermalConformalEpoch out;
  if(!profile() || !(a>=0x1p-512L && a<=1)) { out.status=S::outside_domain; return out; }
  if(!ledger.background()) { out.status=ledger.status; return out; }
  const auto scaled=id.background.scaled_expansion(a);
  out=detail::thermal_conformal_epoch_from_scaled(id.background,scaled,a,k,id.mapping[0].emitted_value);
  if(out.status==S::ok) {
    const auto s=detail::thermal_conformal_no_species_diagnostics(id.background,retained(id),id.mapping,a,out);
    if(s!=S::ok) out.status=s;
  }
  return out;
}
// Nested Simpson age from the actual retained scaled callback, including its
// finite a=0 photon limit. The a=0 opt-in forward/response gate is still absent.
bool age(const FiniteOpacityIdentity &id,W upper,f::Ledger &ledger,W &value,W &difference) {
  W previous=0;
  for(unsigned n=32;n<=65536;n*=2) {
    W sum=0;
    for(unsigned i=0;i<=n;++i) {
      if(!ledger.background()) return false;
      const auto p=id.background.scaled_expansion(upper*i/n);
      if(p.status!=S::ok || !(p.a4_e2>0)) { ledger.status=p.status==S::ok?S::outside_domain:p.status; return false; }
      const W v=detail::thermal_conformal_c_km_s/id.original.model.h0_km_s_mpc/std::sqrt(p.a4_e2);
      if(!std::isfinite(v)) { ledger.status=S::overflow; return false; }
      sum+=(i==0 || i==n?1:i%2?4:2)*v;
    }
    value=upper*sum/(3*n);
    if(n>32) {
      difference=std::abs(value-previous);
      if(difference<=1e-11L*(1+std::abs(value))) return true;
    }
    previous=value;
  }
  ledger.status=S::conditioning_budget_exceeded; return false;
}
// Exact integral of the declared piecewise linear binary64 opacity law.
std::pair<W,W> opacity_and_tau(const SuppliedConformalOpacity &o,W eta) {
  const auto upper=std::upper_bound(o.eta_mpc.begin(),o.eta_mpc.end(),eta);
  const std::size_t j=upper==o.eta_mpc.end()?o.eta_mpc.size()-2:std::size_t(upper-o.eta_mpc.begin()-1);
  const W left=o.eta_mpc[j],right=o.eta_mpc[j+1];
  const W kl=o.differential_opacity_mpc_inverse[j],kr=o.differential_opacity_mpc_inverse[j+1];
  const W k=(kl*(right-eta)+kr*(eta-left))/(right-left);
  W tau=(right-eta)*(k+kr)/2;
  for(std::size_t i=j+1;i+1<o.eta_mpc.size();++i)
    tau+=(W(o.eta_mpc[i+1])-o.eta_mpc[i])*(W(o.differential_opacity_mpc_inverse[i+1])+o.differential_opacity_mpc_inverse[i])/2;
  return {k,tau};
}
bool clock_step(const FiniteOpacityIdentity &id,W k,W h,W a,
                std::array<detail::ThermalConformalEpoch,2> &stages,
                W &next,f::Ledger &ledger,FiniteOpacityAttemptReceipt &receipt) {
  if(!ledger.writes(2))return false;
  std::array<W,2> stage_a{a,a}; W residual=0;
  for(unsigned iteration=0;iteration<24;++iteration) {
    for(unsigned s=0;s<2;++s) {
      stages[s]=epoch(id,stage_a[s],k,ledger);
      if(stages[s].status!=S::ok) { ledger.status=stages[s].status; return false; }
    }
    if(!ledger.writes(4))return false;
    std::array<W,2> revised; residual=0;
    for(unsigned s=0;s<2;++s) {
      revised[s]=a+h*(f::radau_a[s][0]*stage_a[0]*stages[0].hcal+
                      f::radau_a[s][1]*stage_a[1]*stages[1].hcal);
      residual=std::max(residual,std::abs(revised[s]-stage_a[s])/a);
      if(!(revised[s]>=a && revised[s]<=1) || !std::isfinite(revised[s])) { ledger.status=S::outside_domain; return false; }
    }
    if(!ledger.writes(2)) return false;
    stage_a=revised;
    if(residual<=64*std::numeric_limits<W>::epsilon()) {
      for(unsigned s=0;s<2;++s) {
        stages[s]=epoch(id,stage_a[s],k,ledger);
        if(stages[s].status!=S::ok) { ledger.status=stages[s].status; return false; }
      }
      next=stage_a[1];
      receipt.maximum_clock_stage_residual=std::max(receipt.maximum_clock_stage_residual,residual);
      return true;
    }
  }
  receipt.maximum_clock_stage_residual=std::max(receipt.maximum_clock_stage_residual,residual);
  ledger.status=S::conditioning_budget_exceeded; return false;
}
void source_payload(irred::detail::PayloadAccounting &a,const projection::ContinuousCmbSource &s) {
  a.vector(s.k_mpc_inverse);a.vector(s.eta_mpc);a.vector(s.t0);a.vector(s.t1);a.vector(s.t2);a.vector(s.polarization);
  a.string(s.producer_id);a.string(s.signed_mode_id);a.string(s.normalization_id);
}
std::optional<std::size_t> payload(const FiniteOpacitySourceResult &r,const std::vector<double> &grid,
                                 std::size_t scratch) {
  irred::detail::PayloadAccounting a(sizeof(r)+scratch);a.vector(grid);
  if(r.identity) { const auto n=r.identity->retained_payload_bound();if(!n)return {};a.add(*n,1); }
  a.vector(r.attempts);a.vector(r.diagnostics);
  for(const auto &attempt:r.attempts) {
    a.vector(attempt.nodes);a.vector(attempt.final_temperature_tail);a.vector(attempt.final_polarization_tail);
  }
  if(r.source) source_payload(a,*r.source);
  return a.result();
}
bool observe_payload(FiniteOpacitySourceResult &r,const std::vector<double> &grid,std::size_t scratch,
                     const FiniteOpacityPolicy &p) {
  const auto n=payload(r,grid,scratch);
  if(!n) { r.status=S::overflow;return false; }
  if(!r.peak_owned_payload_bound || *n>*r.peak_owned_payload_bound) r.peak_owned_payload_bound=*n;
  if(*n>p.maximum_native_bytes) { r.status=S::work_limit;return false; }
  return true;
}
// Before reserve, charge the old live owners AND the complete requested new
// buffer. Actual capacity is checked immediately afterward; a refused or failed
// acquisition keeps its peak receipt. This excludes allocator metadata/RSS.
template<class T> bool reserve(std::vector<T> &v,std::size_t n,
                              FiniteOpacitySourceResult &r,const std::vector<double> &grid,
                              std::size_t scratch,const FiniteOpacityPolicy &policy) {
  if(n<=v.capacity())return true;
  const auto now=payload(r,grid,scratch);
  if(!now){r.status=S::overflow;return false;}
  std::size_t requested=*now;
  if(!irred::detail::checked_payload_add(requested,n,sizeof(T))){r.status=S::overflow;return false;}
  if(!r.peak_owned_payload_bound || requested>*r.peak_owned_payload_bound)r.peak_owned_payload_bound=requested;
  if(requested>policy.maximum_native_bytes){r.status=S::work_limit;return false;}
  v.reserve(n);
  return observe_payload(r,grid,scratch,policy);
}
bool append_node(const FiniteOpacityIdentity &id,std::size_t ki,std::size_t ti,W eta_mpc,W a,
                 const f::TransportState &state,FiniteOpacityAttemptReceipt &attempt,f::Ledger &ledger) {
  const W k=id.original.k_mpc_inverse[ki];const auto e=epoch(id,a,k,ledger);
  if(e.status!=S::ok) { ledger.status=e.status;return false; }
  const auto [opacity,tau]=opacity_and_tau(id.original.opacity,eta_mpc);
  const W survival=std::exp(-tau),visibility=opacity*survival;
  if(!(survival>0) || !std::isnormal(survival) || !std::isfinite(visibility) ||
     (opacity>0 && (!std::isnormal(visibility) || !(visibility>0)))) { ledger.status=S::overflow;return false; }
  FiniteOpacityNode node;node.eta_index=ti;node.k_index=ki;node.eta_mpc=eta_mpc;node.scale_factor=a;
  node.core=state.core;node.raw_channels=f::raw_source(e,k,opacity,survival,state.core);
  node.psi=f::potential_psi(e,k,state.core);node.phi_prime_mpc_inverse=f::potential_prime(e,k,state.core,node.psi);
  node.hamiltonian_residual=f::hamiltonian(e,k,state.core,node.hamiltonian_term_scale);
  node.normalized_hamiltonian_residual=node.hamiltonian_term_scale>0?
      std::abs(node.hamiltonian_residual)/node.hamiltonian_term_scale:std::abs(node.hamiltonian_residual);
  node.closure_defect=e.closure_defect;node.actual_p_minus_shadow=e.shadow->actual_p_minus_shadow;
  node.conditional_hcal_estimate=e.forward->hcal_absolute_estimate;
  node.survival=survival;node.visibility_mpc_inverse=visibility;
  for(unsigned ch=0;ch<4;++ch) {
    const double emitted=static_cast<double>(node.raw_channels[ch]);
    if(!std::isfinite(emitted) || (node.raw_channels[ch]!=0 && (emitted==0 || !std::isnormal(emitted)))) {
      ledger.status=S::overflow;return false;
    }
    node.channel_cast_loss[ch]=std::abs(node.raw_channels[ch]-W(emitted));
  }
  if(!ledger.writes(11+4+4+12)) return false;
  attempt.nodes.push_back(node);return true;
}
// Slots are acquired once. The mandatory prefix copy after every attempted
// transport solve has reserved write space before that solve is allowed.
void save_tail(const f::TransportState &state,FiniteOpacityAttemptReceipt &attempt,
               unsigned hierarchy,std::size_t ki,f::Ledger &ledger) {
  const std::size_t offset=ki*(hierarchy-2);
  for(unsigned l=3;l<=hierarchy;++l) {
    attempt.final_temperature_tail[offset+l-3]=state.temperature[l];
    attempt.final_polarization_tail[offset+l-3]=state.polarization[l];
  }
  ledger.commit_refusal_writes(2*(hierarchy-2));
}
} // namespace

std::optional<std::size_t> FiniteOpacityIdentity::retained_payload_bound() const noexcept {
  irred::detail::PayloadAccounting a(sizeof(*this));request_payload(a,original);a.vector(seeds);
  a.vector(boundary.photon_monopole);a.vector(boundary.photon_dipole_theta_over_k);
  a.vector(boundary.omitted_temperature_absolute_bound);
  // Empty species are required, but a preexisting empty vector may own capacity.
  a.vector(background.source().species);return a.result();
}
FiniteOpacitySourceProducer::FiniteOpacitySourceProducer(FiniteOpacitySourceProducer &&o) noexcept
    :status_(std::exchange(o.status_,S::invalid_input)),identity_(std::move(o.identity_)) {}
FiniteOpacitySourceProducer &FiniteOpacitySourceProducer::operator=(FiniteOpacitySourceProducer &&o) noexcept {
  if(this!=&o) {status_=std::exchange(o.status_,S::invalid_input);identity_=std::move(o.identity_);}return *this;
}
std::optional<std::size_t> FiniteOpacitySourceProducer::retained_payload_bound() const noexcept {
  irred::detail::PayloadAccounting a(sizeof(*this));
  if(identity_) {const auto n=identity_->retained_payload_bound();if(!n)return {};a.add(*n,1);}return a.result();
}
FiniteOpacitySourceProducer prepare_finite_opacity_source(FiniteOpacityRequest &&request,FiniteOpacityPolicy policy) {
  FiniteOpacitySourceProducer out;out.status_=request_valid(request,policy);
  if(out.status_!=S::ok)return out;
  try {
    auto id=std::make_shared<FiniteOpacityIdentity>();id->original=std::move(request);id->preparation_policy=policy;
    id->model_id=finite_opacity_model_id;id->mode_id=finite_opacity_mode_id;
    id->method_id=finite_opacity_method_id;id->arithmetic_id=finite_opacity_arithmetic_id;
    // Retain the acquired original even on mapper/clock/seed refusal.
    out.identity_=id;
    const auto map=map_thermal_physical_model(id->original.model);
    out.status_=map.status;if(!map.model || !map.scalar_witnesses)return out;
    id->mapping=*map.scalar_witnesses;
    id->background=prepare_thermal_background(*map.model);
    out.status_=id->background.status();if(out.status_!=S::ok)return out;
    const auto receipt=detail::ThermalRetainedCoefficientAccess::capture(id->background);
    if(!receipt){out.status_=S::outside_domain;return out;}
    id->lambda_retained=receipt->lambda_retained;id->lambda_emitted=receipt->lambda_emitted;
    id->lambda_getter_signed_loss=receipt->lambda_getter_signed_loss;
    id->critical_density_ev4=receipt->critical_density_ev4;
    id->species_normalization_error=receipt->species_normalization_error;
    f::Ledger ledger{id->preparation_work,policy};
    if(!age(*id,id->original.initial_scale_factor,ledger,id->clock.initial_age_mpc,id->clock.initial_quadrature_difference_mpc) ||
       !age(*id,id->original.observer_scale_factor,ledger,id->clock.observer_age_mpc,id->clock.observer_quadrature_difference_mpc)) {
      out.status_=ledger.status;return out;
    }
    id->clock.initial_eta_mismatch_mpc=W(id->original.opacity.eta_mpc.front())-id->clock.initial_age_mpc;
    id->clock.observer_eta_mismatch_mpc=W(id->original.opacity.eta_mpc.back())-id->clock.observer_age_mpc;
    const auto initial=epoch(*id,id->original.initial_scale_factor,id->original.k_mpc_inverse.front(),ledger);
    const auto observer=epoch(*id,id->original.observer_scale_factor,id->original.k_mpc_inverse.front(),ledger);
    if(initial.status!=S::ok || observer.status!=S::ok) {out.status_=initial.status!=S::ok?initial.status:observer.status;return out;}
    // Conditional endpoint consistency only. These tolerances do not certify
    // the integrated clock response or assign its absent source error to zero.
    const W initial_limit=policy.state_absolute_tolerance/(5*initial.hcal)+id->clock.initial_quadrature_difference_mpc;
    const W observer_limit=policy.state_absolute_tolerance/(5*observer.hcal)+id->clock.observer_quadrature_difference_mpc;
    id->clock.endpoint_clock_consistent=std::abs(id->clock.initial_eta_mismatch_mpc)<=initial_limit &&
        std::abs(id->clock.observer_eta_mismatch_mpc)<=observer_limit;
    if(!id->clock.endpoint_clock_consistent){out.status_=S::conditioning_budget_exceeded;return out;}
    id->seeds.reserve(id->original.k_mpc_inverse.size());
    for(double k:id->original.k_mpc_inverse) {
      if(!ledger.writes(33)) {out.status_=ledger.status;return out;}
      auto e=initial;e.x2=(W(k)/e.hcal)*(W(k)/e.hcal);
      const auto seed=f::finite_seed(e,k);FiniteOpacitySeed saved;f::Core emitted{};
      for(unsigned i=0;i<11;++i) {
        saved.emitted_core[i]=static_cast<double>(seed[i]);
        if(!std::isfinite(saved.emitted_core[i]) || (seed[i]!=0 && (saved.emitted_core[i]==0 || !std::isnormal(saved.emitted_core[i])))) {
          out.status_=S::overflow;return out;
        }
        emitted[i]=saved.emitted_core[i];saved.absolute_cast_loss[i]=std::abs(seed[i]-emitted[i]);
      }
      W scale=0;saved.hamiltonian_residual=f::hamiltonian(e,k,emitted,scale);
      saved.zeta_residual=-emitted[f::phi]+emitted[f::db]/3-1;
      id->seeds.push_back(saved);
    }
    const auto op=opacity_and_tau(id->original.opacity,id->original.opacity.eta_mpc.front());
    auto &boundary=id->boundary;boundary.eta_i_mpc=id->original.opacity.eta_mpc.front();
    boundary.tau_i=op.second;boundary.survival_i=std::exp(-op.second);
    boundary.photon_monopole.reserve(id->seeds.size());boundary.photon_dipole_theta_over_k.reserve(id->seeds.size());
    boundary.omitted_temperature_absolute_bound.reserve(id->seeds.size());
    for(const auto &seed:id->seeds) {
      if(!ledger.writes(3)){out.status_=ledger.status;return out;}
      const W m=seed.emitted_core[f::dg]/4,u=seed.emitted_core[f::ug];
      boundary.photon_monopole.push_back(m);boundary.photon_dipole_theta_over_k.push_back(u);
      boundary.omitted_temperature_absolute_bound.push_back(boundary.survival_i*(std::abs(m)+std::abs(u)));
    }
    const auto bytes=out.retained_payload_bound();
    if(!bytes)out.status_=S::overflow;else if(*bytes>policy.maximum_native_bytes)out.status_=S::work_limit;
  } catch(const std::bad_alloc &) {out.status_=S::work_limit;}
  catch(const std::length_error &) {out.status_=S::overflow;}
  return out;
}

FiniteOpacitySourceResult FiniteOpacitySourceProducer::produce(FiniteOpacityPolicy policy) const {
  FiniteOpacitySourceResult out;out.identity=identity_;out.status=status_;out.policy=policy;
  if(!identity_)return out;
  if(const auto bytes=identity_->retained_payload_bound()) {
    std::size_t retained_bytes=*bytes;
    if(irred::detail::checked_payload_add(retained_bytes,1,sizeof(out)))out.peak_owned_payload_bound=retained_bytes;
  }
  out.work=identity_->preparation_work;
  if(!profile()) {out.status=S::outside_domain;return out;}
  if(status_!=S::ok)return out;
  out.status=request_valid(identity_->original,policy);if(out.status!=S::ok)return out;
  if(out.work.background_clock_calls>policy.maximum_background_clock_calls ||
     out.work.destination_writes>policy.maximum_destination_writes) {out.status=S::work_limit;return out;}
  f::Ledger ledger{out.work,policy};std::vector<double> grid;
  constexpr std::size_t scratch_bytes=sizeof(f::TransportScratch)+sizeof(f::TransportState)+
      4*sizeof(detail::ThermalConformalEpoch)+4*sizeof(f::Core)+
      2*sizeof(detail::ThermalRetainedCoefficientWitness)+sizeof(f::Ledger)+4096;
  try {
    const auto &id=*identity_;const auto &o=id.original.opacity;const std::size_t nk=id.original.k_mpc_inverse.size();
    if(!reserve(grid,4*(o.eta_mpc.size()-1)+1,out,grid,scratch_bytes,policy))return out;
    if(!ledger.writes(4*(o.eta_mpc.size()-1)+1)){out.status=ledger.status;return out;}
    for(std::size_t i=0;i+1<o.eta_mpc.size();++i) for(unsigned q=0;q<4;++q)
      grid.push_back(static_cast<double>(W(o.eta_mpc[i])+(W(o.eta_mpc[i+1])-o.eta_mpc[i])*q/4));
    grid.push_back(o.eta_mpc.back());
    for(std::size_t i=1;i<grid.size();++i)if(!(grid[i]>grid[i-1])){out.status=S::overflow;return out;}
    // Nine trajectories on the SAME finest producer grid. Every original
    // opacity knot is preserved; its output interpolation error remains absent.
    if(!reserve(out.attempts,9,out,grid,scratch_bytes,policy))return out;
    if(!observe_payload(out,grid,scratch_bytes,policy))return out;
    for(unsigned lmax:{48u,96u,192u}) for(unsigned refinement=0;refinement<3;++refinement) {
      out.attempts.emplace_back();auto &attempt=out.attempts.back();
      attempt.hierarchy=lmax;attempt.time_refinement=refinement;attempt.status=S::ok;
      if(!reserve(attempt.nodes,grid.size()*nk,out,grid,scratch_bytes,policy) ||
         !reserve(attempt.final_temperature_tail,(lmax-2)*nk,out,grid,scratch_bytes,policy) ||
         !reserve(attempt.final_polarization_tail,(lmax-2)*nk,out,grid,scratch_bytes,policy)){
        attempt.status=out.status;return out;
      }
      // Initialize all endpoint slots before trajectory work. Mandatory copies
      // reserve separate write space before each full transport solve.
      if(!ledger.writes(2*(lmax-2)*nk)){
        attempt.status=out.status=ledger.status;return out;
      }
      attempt.final_temperature_tail.resize((lmax-2)*nk);
      attempt.final_polarization_tail.resize((lmax-2)*nk);
      if(!observe_payload(out,grid,scratch_bytes,policy)) {attempt.status=out.status;return out;}
      for(std::size_t ki=0;ki<nk;++ki) {
        const W k=id.original.k_mpc_inverse[ki];W a=id.original.initial_scale_factor,eta_mpc=grid.front();
        // Both full fixed-capacity model buffers are initialized once per k;
        // unused hierarchy slots still belong to the actual scratch owner.
        constexpr std::size_t scratch_initializations=2*193*10+22*22+44+2*(11+2*193);
        if(!ledger.writes(scratch_initializations+11+2*193+11)){attempt.status=out.status=ledger.status;return out;}
        f::TransportState state;f::TransportScratch scratch;
        attempt.reached_wavenumbers=ki+1;
        for(unsigned j=0;j<11;++j)state.core[j]=id.seeds[ki].emitted_core[j];
        if(!append_node(id,ki,0,eta_mpc,a,state,attempt,ledger)) {attempt.status=out.status=ledger.status;return out;}
        for(std::size_t ti=1;ti<grid.size();++ti) {
          while(eta_mpc<grid[ti]) {
            ++attempt.attempted_steps;
            if(!ledger.step()){attempt.status=out.status=ledger.status;return out;}
            const auto current=epoch(id,a,k,ledger);
            if(current.status!=S::ok){attempt.status=out.status=current.status;return out;}
            const W factor=std::ldexp(1.L,-int(refinement));
            const W h=std::min({W(grid[ti])-eta_mpc,W(policy.maximum_eta_step_mpc)*factor,
                              .1L*factor/k,W(policy.maximum_log_a_step)*factor/current.hcal});
            if(!(h>0) || eta_mpc+h==eta_mpc) {attempt.status=out.status=S::overflow;return out;}
            std::array<detail::ThermalConformalEpoch,2> stages;W next=0;
            if(!clock_step(id,k,h,a,stages,next,ledger,attempt)) {attempt.status=out.status=ledger.status;return out;}
            const f::Pair opacity{opacity_and_tau(o,eta_mpc+h/3).first,opacity_and_tau(o,eta_mpc+h).first};
            if(!ledger.reserve_refusal_writes(2*(lmax-2))){attempt.status=out.status=ledger.status;return out;}
            const bool solved=f::radau_step(stages,opacity,lmax,eta_mpc,h,k,state,scratch,ledger,attempt);
            save_tail(state,attempt,lmax,ki,ledger);
            if(!solved) {
              attempt.status=out.status=ledger.status;return out;
            }
            a=next;eta_mpc=std::min(W(grid[ti]),eta_mpc+h);++attempt.completed_steps;
          }
          if(!append_node(id,ki,ti,eta_mpc,a,state,attempt,ledger)) {attempt.status=out.status=ledger.status;return out;}
        }
        attempt.maximum_observer_scale_factor_mismatch=std::max(attempt.maximum_observer_scale_factor_mismatch,
            std::abs(a-W(id.original.observer_scale_factor)));
      }
      attempt.observer_scale_factor_consistent=attempt.maximum_observer_scale_factor_mismatch<=
          (policy.state_absolute_tolerance+policy.relative_tolerance*id.original.observer_scale_factor)/5;
    }
    // A full computed raw source is useful before integrated error/reference
    // closure. Complete nine-attempt prefixes are required even for that export.
    const auto &finest=out.attempts[8];out.source.emplace();auto &source=*out.source;
    if(!reserve(source.k_mpc_inverse,nk,out,grid,scratch_bytes,policy) ||
       !reserve(source.eta_mpc,grid.size(),out,grid,scratch_bytes,policy) ||
       !reserve(source.t0,grid.size()*nk,out,grid,scratch_bytes,policy) ||
       !reserve(source.t1,grid.size()*nk,out,grid,scratch_bytes,policy) ||
       !reserve(source.t2,grid.size()*nk,out,grid,scratch_bytes,policy) ||
       !reserve(source.polarization,grid.size()*nk,out,grid,scratch_bytes,policy)) {out.source.reset();return out;}
    if(!ledger.writes(nk+4*grid.size()*nk)){out.status=ledger.status;out.source.reset();return out;}
    source.k_mpc_inverse.assign(id.original.k_mpc_inverse.begin(),id.original.k_mpc_inverse.end());
    source.t0.resize(grid.size()*nk);source.t1.resize(grid.size()*nk);
    source.t2.resize(grid.size()*nk);source.polarization.resize(grid.size()*nk);
    source.observer_eta_mpc=o.eta_mpc.back();source.producer_id=std::string(finite_opacity_model_id);
    source.signed_mode_id=std::string(finite_opacity_mode_id);source.normalization_id="finite-local-zeta-one/computed-unadmitted/v1";
    if(!reserve(out.diagnostics,grid.size()*nk*4,out,grid,scratch_bytes,policy)){out.source.reset();return out;}
    if(!observe_payload(out,grid,scratch_bytes,policy)){out.source.reset();return out;}
    for(W eta:grid) {
      const double emitted=static_cast<double>(eta);
      if(!std::isfinite(emitted) || (!source.eta_mpc.empty() && !(emitted>source.eta_mpc.back()))) {
        out.status=S::overflow;out.source.reset();return out;
      }
      if(!ledger.writes(1)){out.status=ledger.status;out.source.reset();return out;}
      source.eta_mpc.push_back(emitted);
    }
    for(std::size_t ki=0;ki<nk;++ki) for(std::size_t ti=0;ti<grid.size();++ti) {
      const auto index=ki*grid.size()+ti;const auto &node=finest.nodes[index];
      for(unsigned ch=0;ch<4;++ch) {
        const W value=node.raw_channels[ch];FiniteOpacityChannelDiagnostic d;
        d.eta_index=ti;d.k_index=ki;d.channel=ch;
        d.allocation=policy.channel_absolute_tolerance_mpc_inverse+policy.relative_tolerance*std::abs(value);
        d.time_difference=std::abs(value-out.attempts[7].nodes[index].raw_channels[ch]);
        d.hierarchy_difference=std::abs(value-out.attempts[5].nodes[index].raw_channels[ch]);
        d.previous_time_difference=std::abs(out.attempts[7].nodes[index].raw_channels[ch]-out.attempts[6].nodes[index].raw_channels[ch]);
        d.previous_hierarchy_difference=std::abs(out.attempts[5].nodes[index].raw_channels[ch]-out.attempts[2].nodes[index].raw_channels[ch]);
        d.measured_cast_loss=node.channel_cast_loss[ch];out.diagnostics.push_back(d);
        if(!ledger.writes(1)){out.status=ledger.status;out.source.reset();return out;}
        auto &v=ch==0?source.t0:ch==1?source.t1:ch==2?source.t2:source.polarization;
        v[ti*nk+ki]=static_cast<double>(value);
      }
    }
    out.status=S::conditioning_budget_exceeded; // no missing error assigned zero
    out.source_numerically_admitted=false;
    observe_payload(out,grid,0,policy);
  } catch(const std::bad_alloc &) {out.status=S::work_limit;observe_payload(out,grid,scratch_bytes,policy);out.source.reset();}
  catch(const std::length_error &) {out.status=S::overflow;observe_payload(out,grid,scratch_bytes,policy);out.source.reset();}
  return out;
}
} // namespace irred::cosmology
