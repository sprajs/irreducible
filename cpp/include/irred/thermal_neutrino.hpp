#pragma once
#include "irred/numerics.hpp"
#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>
namespace irred::cosmology {
// Collisionless f(q)=1/(exp(q)+1), q=p*a/Tnu0, in natural hbar=c=kB=1
// units. Statistical weight counts populated states explicitly: e.g. g=2
// for one neutrino plus antineutrino. No temperature/Neff/mass-sum mapping.
struct ThermalSpecies {
  double mass_ev, temperature_today_ev, statistical_weight;
};
enum class ThermalMomentumMethod { direct_adaptive, nested_clenshaw_curtis };
struct ThermalPolicy {
  // Applies separately to I_rho/hypot(1,y) and I_pressure*hypot(1,y).
  double absolute_tolerance = 1e-12, relative_tolerance = 2e-12;
  std::size_t maximum_callbacks_per_evaluation = 200000;
  std::size_t maximum_total_callbacks = 4000000;
  unsigned maximum_depth = 30;
  std::size_t maximum_points = 4096, maximum_species = 16;
  std::size_t maximum_native_bytes = 16 * 1024 * 1024;
  // Preparation and later retained-owner evaluations must select the same method.
  ThermalMomentumMethod momentum_method = ThermalMomentumMethod::direct_adaptive;
};
struct ThermalMoments {
  numerics::Status status = numerics::Status::invalid_input;
  // I_rho=integral q^2 sqrt(q^2+y^2) f dq;
  // I_pressure=integral q^4/[3 sqrt(q^2+y^2)] f dq.
  std::optional<double> rho_moment, pressure_moment;
  double rho_error_estimate = 0, pressure_error_estimate = 0;
  std::size_t callbacks = 0;
};
ThermalMoments evaluate_thermal_moments(double y, ThermalPolicy = {});
struct ThermalSpeciesState {
  numerics::Status status = numerics::Status::invalid_input;
  // Natural-unit energy density and pressure, eV^4 (not eV/m^3).
  std::optional<double> rho_ev4, pressure_ev4;
  double rho_error_estimate_ev4 = 0, pressure_error_estimate_ev4 = 0;
  std::size_t callbacks = 0;
};
ThermalSpeciesState evaluate_thermal_species(const ThermalSpecies &, double a,
                                             ThermalPolicy = {});
struct ThermalFlatModel {
  double h0_km_s_mpc, omega_gamma, omega_massless_nonphoton, omega_b, omega_cdm;
  std::vector<ThermalSpecies> species;
};
// Source physical densities are omega_i = Omega_i h^2, h=H0/100.
// Photon energy density derives from the explicitly supplied Kelvin
// temperature.
struct ThermalPhysicalSpecies {
  double mass_ev, temperature_today_kelvin, statistical_weight;
};
struct ThermalPhysicalModel {
  double h0_km_s_mpc, physical_baryon_density, physical_cdm_density;
  double tcmb_kelvin, physical_massless_nonphoton_density;
  std::vector<ThermalPhysicalSpecies> species;
};
// Once-captured deterministic arithmetic witnesses of the original map.
// Scalar order is photon, baryon, CDM, other massless density. These are not
// physical-input uncertainties; each consuming model owns propagation.
struct ThermalScalarMapWitness {
  long double wide_value = 0, wide_operation_estimate = 0;
  long double measured_absolute_cast_loss = 0;
  double emitted_value = 0;
};
struct ThermalPhysicalMapping {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<ThermalFlatModel> model;
  std::optional<std::array<ThermalScalarMapWitness, 4>> scalar_witnesses;
};
ThermalPhysicalMapping map_thermal_physical_model(const ThermalPhysicalModel &,
                                                  ThermalPolicy = {});
// Wide internal numerical coordinate, P(a)=a^4 E(a)^2. The finite a=0
// radiation limit is meaningful; an exactly radiation-free source has P(0)=0.
// E(0) is not evaluated or fabricated. Nonzero wide coordinates/diagnostics
// must remain normal; positive-a underflow cannot masquerade as exact zero.
struct ThermalScaledExpansion {
  numerics::Status status = numerics::Status::invalid_input;
  long double a4_e2 = 0, error_estimate = 0;
  std::size_t callbacks = 0;
};
struct ThermalBackgroundValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double error_estimate = 0;
};
inline constexpr unsigned thermal_e = 1, thermal_h = 2;
struct ThermalBackgroundRow {
  double scale_factor = 0;
  ThermalBackgroundValue e, h_km_s_mpc;
  std::size_t callbacks = 0;
};
struct ThermalBackgroundBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested_outputs = 0;
  std::vector<ThermalBackgroundRow> rows;
  std::size_t callbacks = 0;
};
class ThermalBackground {
public:
  ThermalBackground() = default;
  ThermalBackground(const ThermalBackground &) = default;
  ThermalBackground &operator=(const ThermalBackground &) = default;
  // Moving transfers the complete physical owner and invalidates the source.
  // Self move preserves the current owner; copies retain independent storage.
  ThermalBackground(ThermalBackground &&) noexcept;
  ThermalBackground &operator=(ThermalBackground &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const ThermalFlatModel &source() const noexcept { return source_; }
  ThermalMomentumMethod momentum_method() const noexcept { return method_; }
  std::optional<double> omega_species_today() const noexcept;
  std::optional<double> omega_lambda() const noexcept;
  std::size_t preparation_callbacks() const noexcept { return callbacks_; }
  ThermalScaledExpansion scaled_expansion(long double a,
                                          ThermalPolicy = {}) const;
  ThermalBackgroundBatch evaluate(std::span<const double> scale_factors,
                                  unsigned requested_outputs,
                                  ThermalPolicy = {}) const;

private:
  numerics::Status status_ = numerics::Status::invalid_input;
  ThermalFlatModel source_{};
  ThermalMomentumMethod method_ = ThermalMomentumMethod::direct_adaptive;
  long double critical_ev4_ = 0, omega_species_ = 0, lambda_ = 0;
  long double normalization_error_ = 0;
  std::size_t callbacks_ = 0;
  friend ThermalBackground prepare_thermal_background(const ThermalFlatModel &,
                                                      ThermalPolicy);
};
ThermalBackground prepare_thermal_background(const ThermalFlatModel &,
                                             ThermalPolicy = {});
// Conservative simultaneous owned payload, excluding input owners, allocator
// metadata, stack recursion and RSS. Hard domain: <=16 explicit species.
std::optional<std::size_t>
thermal_background_payload_bound(std::size_t points,
                                 std::size_t species) noexcept;
inline constexpr std::string_view thermal_neutrino_model_id =
    "flat-collisionless-zero-chemical-potential-thermal-FD-relic-lambda/v1";
inline constexpr std::string_view thermal_neutrino_method_id =
    "scaled-adaptive-direct-momentum-exponential-tail/v1";
inline constexpr std::string_view thermal_neutrino_arithmetic_id =
    "thermal-FD/binary64-quadrature-wide-scaling/v1";
inline constexpr std::string_view thermal_neutrino_constants_id =
    "SI2019-exact-h-c-kB-eV-IAU2012-AU-CODATA2018-G-fixed";
} // namespace irred::cosmology

namespace irred::cosmology {
// The opt-in method includes bounded refinement and charged direct fallback.
constexpr std::string_view thermal_momentum_method_id(ThermalMomentumMethod method) {
  switch (method) {
  case ThermalMomentumMethod::direct_adaptive: return thermal_neutrino_method_id;
  case ThermalMomentumMethod::nested_clenshaw_curtis:
    return "scaled-nested-clenshaw-curtis-direct-fallback-exponential-tail/v1";
  }
  return {};
}
constexpr std::string_view thermal_momentum_arithmetic_id(ThermalMomentumMethod method) {
  switch (method) {
  case ThermalMomentumMethod::direct_adaptive: return thermal_neutrino_arithmetic_id;
  case ThermalMomentumMethod::nested_clenshaw_curtis:
    return "thermal-FD/binary64-shared-nodes-wide-sums/v1";
  }
  return {};
}
} // namespace irred::cosmology
