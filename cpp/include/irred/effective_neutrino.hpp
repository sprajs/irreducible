#pragma once
#include "irred/thermal_observables.hpp"
#include <string_view>
namespace irred::cosmology {
// Degeneracy counts neutrino/antineutrino pairs: native FD g=2*degeneracy.
// Every value is explicit. No hierarchy, mass-sum or 93.14-eV conversion.
struct EffectiveNeutrinoSpecies {
  double mass_ev, temperature_ratio, degeneracy;
};
struct EffectiveNeutrinoModel {
  double h0_km_s_mpc, physical_baryon_density, physical_cdm_density;
  double tcmb_kelvin, effective_relativistic_species;
  std::vector<EffectiveNeutrinoSpecies> species;
};
struct EffectiveNeutrinoDiagnostics {
  long double explicit_species_early_neff = 0, remaining_massless_neff = 0;
  long double partition_error_estimate = 0;
  double realized_early_neff = 0;
  long double realized_early_neff_error_estimate = 0;
  std::optional<double> physical_relic_density_today, omega_lambda;
  // Actual physical-map scalar witnesses in photon,baryon,CDM,ur order.
  std::array<ThermalScalarMapWitness,4> scalar_witnesses{};
  ThermalStressSourceEstimate stress_source_estimate;
};
class EffectiveNeutrinoState {
public:
  EffectiveNeutrinoState() = default;
  EffectiveNeutrinoState(const EffectiveNeutrinoState &) = default;
  EffectiveNeutrinoState &operator=(const EffectiveNeutrinoState &) = default;
  EffectiveNeutrinoState(EffectiveNeutrinoState &&) noexcept;
  EffectiveNeutrinoState &operator=(EffectiveNeutrinoState &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const EffectiveNeutrinoModel &source() const noexcept { return source_; }
  const ThermalPhysicalModel &physical_source() const noexcept { return physical_; }
  const ThermalBackground &background() const noexcept { return background_; }
  const EffectiveNeutrinoDiagnostics &diagnostics() const noexcept { return diagnostics_; }
  ThermalStressBatch evaluate_stress_energy(std::span<const double>,
                                           ThermalStressPolicy = {}) const;
  // The returned request composes existing distance/ruler and BAO consumers.
  // It retains supplied-drag identity; this state predicts no drag epoch.
  std::optional<ThermalObservableRequest> observable_request(
      double z_drag, std::string drag_origin, std::string source_origin) const;
private:
  numerics::Status status_ = numerics::Status::invalid_input;
  EffectiveNeutrinoModel source_{};
  ThermalPhysicalModel physical_{};
  ThermalBackground background_;
  EffectiveNeutrinoDiagnostics diagnostics_{};
  friend EffectiveNeutrinoState prepare_effective_neutrino_state(
      const EffectiveNeutrinoModel &, ThermalPolicy);
};
EffectiveNeutrinoState prepare_effective_neutrino_state(
    const EffectiveNeutrinoModel &, ThermalPolicy = {});
std::optional<std::size_t> effective_neutrino_payload_bound(std::size_t species) noexcept;
inline constexpr std::string_view effective_neutrino_model_id =
    "flat-effective-thermal-FD-total-Neff-partition-lambda/v1";
inline constexpr std::string_view effective_neutrino_mapping_id =
    "total-Neff-minus-explicit-FD-pairs-temperature-fourth-power-SI/v1";
} // namespace irred::cosmology
