#include "irred/effective_neutrino.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <utility>
namespace irred::cosmology {
namespace {
using S=numerics::Status;
constexpr long double eps=std::numeric_limits<double>::epsilon();
bool representable(long double x) {
  return std::isfinite(x) && (x==0 ||
      (x>=std::numeric_limits<double>::min() && x<=std::numeric_limits<double>::max()));
}
}
std::optional<std::size_t> effective_neutrino_payload_bound(std::size_t count) noexcept {
  irred::detail::PayloadAccounting bytes(sizeof(EffectiveNeutrinoState)+sizeof(ThermalPhysicalMapping));
  bytes.add(count,2*sizeof(EffectiveNeutrinoSpecies)+3*sizeof(ThermalPhysicalSpecies)+3*sizeof(ThermalSpecies));
  const auto nested=thermal_background_payload_bound(0,count);
  if (!nested) return {};
  bytes.add(1,*nested);
  return bytes.result();
}
EffectiveNeutrinoState prepare_effective_neutrino_state(
    const EffectiveNeutrinoModel &model, ThermalPolicy policy) {
  EffectiveNeutrinoState out;
  if (std::fegetround()!=FE_TONEAREST || std::numeric_limits<long double>::digits<64 ||
      model.species.size()>16) return out;
  const auto bytes=effective_neutrino_payload_bound(model.species.size());
  if (!bytes) return out;
  if (model.species.size()>policy.maximum_species || *bytes>policy.maximum_native_bytes) {
    out.status_=S::work_limit; return out;
  }
  for (double x:{model.h0_km_s_mpc,model.physical_baryon_density,model.physical_cdm_density,
                 model.tcmb_kelvin,model.effective_relativistic_species})
    if (!std::isfinite(x)) { out.status_=S::nonfinite_input; return out; }
  if (!(model.h0_km_s_mpc>0) || !(model.tcmb_kelvin>0) || model.physical_baryon_density<0 ||
      model.physical_cdm_density<0 || model.effective_relativistic_species<0) {
    out.status_=S::outside_domain; return out;
  }
  // Reuse the authoritative photon conversion; no copied blackbody equation.
  ThermalPhysicalModel physical{model.h0_km_s_mpc,model.physical_baryon_density,
      model.physical_cdm_density,model.tcmb_kelvin,0,{}};
  const auto photon_map=map_thermal_physical_model(physical,policy);
  if (photon_map.status!=S::ok) { out.status_=photon_map.status; return out; }
  const long double r0=std::cbrt(4.L/11), neutrino_photon_ratio=7.L/8*std::pow(r0,4);
  long double explicit_neff=0;
  physical.species.reserve(model.species.size());
  for (const auto &species:model.species) {
    for (double x:{species.mass_ev,species.temperature_ratio,species.degeneracy})
      if (!std::isfinite(x)) { out.status_=S::nonfinite_input; return out; }
    if (species.mass_ev<0 || !(species.temperature_ratio>0) || !(species.degeneracy>0)) {
      out.status_=S::outside_domain; return out;
    }
    const long double n=species.degeneracy*std::pow(static_cast<long double>(species.temperature_ratio)/r0,4);
    const long double temperature=static_cast<long double>(species.temperature_ratio)*model.tcmb_kelvin;
    const long double weight=2.L*species.degeneracy;
    if (!representable(n) || !representable(temperature) || !representable(weight)) {
      out.status_=S::outside_domain; return out;
    }
    explicit_neff+=n;
    physical.species.push_back({species.mass_ev,static_cast<double>(temperature),static_cast<double>(weight)});
  }
  if (!std::isfinite(explicit_neff) || explicit_neff>model.effective_relativistic_species) {
    out.status_=S::outside_domain; return out;
  }
  auto &diagnostic=out.diagnostics_;
  diagnostic.explicit_species_early_neff=explicit_neff;
  diagnostic.remaining_massless_neff=static_cast<long double>(model.effective_relativistic_species)-explicit_neff;
  diagnostic.partition_error_estimate=128*std::numeric_limits<long double>::epsilon()*
      (model.effective_relativistic_species+explicit_neff);
  const long double h=static_cast<long double>(model.h0_km_s_mpc)/100, h2=h*h;
  const auto &photon_witness=(*photon_map.scalar_witnesses)[0];
  const long double omega_ur=photon_witness.wide_value*h2*neutrino_photon_ratio*
      diagnostic.remaining_massless_neff;
  if (!representable(omega_ur)) { out.status_=S::outside_domain; return out; }
  physical.physical_massless_nonphoton_density=static_cast<double>(omega_ur);
  const auto mapping=map_thermal_physical_model(physical,policy);
  if (mapping.status!=S::ok) { out.status_=mapping.status; return out; }
  // Acquire source owners before quadrature so failed preparation work survives.
  out.source_=model; out.physical_=std::move(physical);
  diagnostic.scalar_witnesses=*mapping.scalar_witnesses;
  diagnostic.stress_source_estimate.relative_component_estimate=128*eps;
  diagnostic.stress_source_estimate.radiation_absolute_estimate=
      photon_witness.wide_value*neutrino_photon_ratio*diagnostic.partition_error_estimate+
      std::abs(omega_ur-static_cast<long double>(out.physical_.physical_massless_nonphoton_density))/h2;
  out.background_=prepare_thermal_background(*mapping.model,policy);
  out.status_=out.background_.status();
  if (out.status_!=S::ok) return out;
  // Measure the emitted early radiation coordinate from the actual shared
  // retained owner (including the final eV casts), not the ideal source formula.
  const auto early=out.background_.scaled_expansion(0,policy);
  if (early.status!=S::ok) { out.status_=early.status; return out; }
  const long double photon=mapping.model->omega_gamma;
  const long double realized=(early.a4_e2-photon)/(photon*neutrino_photon_ratio);
  if (!representable(realized) || (model.effective_relativistic_species>0 && realized==0)) {
    out.status_=S::outside_domain; return out;
  }
  diagnostic.realized_early_neff=static_cast<double>(realized);
  diagnostic.realized_early_neff_error_estimate=early.error_estimate/(photon*neutrino_photon_ratio)+
      128*eps*(1+std::abs(realized))+diagnostic.partition_error_estimate+
      std::abs(realized-static_cast<long double>(diagnostic.realized_early_neff));
  const auto species_today=out.background_.omega_species_today();
  if (!species_today) { out.status_=S::outside_domain; return out; }
  const long double relic_omega=*species_today*h2;
  if (!representable(relic_omega)) { out.status_=S::outside_domain; return out; }
  diagnostic.physical_relic_density_today=static_cast<double>(relic_omega);
  diagnostic.omega_lambda=out.background_.omega_lambda();
  return out;
}
EffectiveNeutrinoState::EffectiveNeutrinoState(EffectiveNeutrinoState &&other) noexcept {
  *this=std::move(other);
}
EffectiveNeutrinoState &EffectiveNeutrinoState::operator=(EffectiveNeutrinoState &&other) noexcept {
  if (this==&other) return *this;
  status_=other.status_; source_=std::move(other.source_); physical_=std::move(other.physical_);
  background_=std::move(other.background_); diagnostics_=other.diagnostics_;
  other.status_=S::invalid_input; other.source_={}; other.physical_={}; other.diagnostics_={};
  return *this;
}
ThermalStressBatch EffectiveNeutrinoState::evaluate_stress_energy(
    std::span<const double> factors,ThermalStressPolicy policy) const {
  if (status_!=S::ok) { ThermalStressBatch out; out.status=status_; return out; }
  return background_.evaluate_stress_energy(factors,policy,diagnostics_.stress_source_estimate);
}
std::optional<ThermalObservableRequest> EffectiveNeutrinoState::observable_request(
    double z_drag,std::string drag_origin,std::string source_origin) const {
  if (status_!=S::ok || !std::isfinite(z_drag) || z_drag<0 || drag_origin.empty() || source_origin.empty()) return {};
  return ThermalObservableRequest{physical_,z_drag,std::move(drag_origin),std::move(source_origin)};
}
} // namespace irred::cosmology
