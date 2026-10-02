#pragma once
#include "irred/numerics.hpp"
#include <string_view>
namespace irred::atomic {
// Fixed central atomic assets; their physical uncertainties are excluded from
// numerical diagnostics. No external implementation or partition table.
inline constexpr long double hydrogen_electron_mass_kg = 9.1093837139e-31L;
inline constexpr long double hydrogen_electron_mass_uncertainty_kg = 2.8e-40L;
inline constexpr long double hydrogen_ionization_energy_ev = 13.598434599702L;
inline constexpr long double hydrogen_ionization_energy_uncertainty_ev =
    1.2e-11L;
inline constexpr std::string_view hydrogen_equilibrium_model_id =
    "homogeneous_ground_state_pure_hydrogen_saha_supplied_temperature_density";
inline constexpr std::string_view hydrogen_equilibrium_assets_id =
    "SI2019-exact-h-kB-eV-CODATA2022-me-NIST-ASD5.12-H1-fixed";
inline constexpr std::string_view hydrogen_equilibrium_method_id =
    "log_saha_separate_stable_positive_quadratic_fractions";
inline constexpr std::string_view hydrogen_equilibrium_arithmetic_id =
    "longdouble64-cpu-nearest-strict-fp-v1";
enum HydrogenOutput : unsigned { ionized_fraction = 1, neutral_fraction = 2 };
struct HydrogenState {
  double temperature_kelvin = 0;
  // Physical total hydrogen nuclei density: protons + neutral atoms, m^-3.
  double hydrogen_nuclei_per_cubic_metre = 0;
};
struct HydrogenFraction {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  // Empirical libm/arithmetic/cast diagnostic; never a universal proven bound.
  double relative_arithmetic_estimate = 0;
};
struct HydrogenRow {
  HydrogenState source;
  numerics::Status admission_status = numerics::Status::invalid_input;
  std::size_t solves = 0;
  HydrogenFraction ionized, neutral;
};
struct HydrogenPolicy {
  unsigned requested_outputs = ionized_fraction | neutral_fraction;
  std::size_t maximum_rows = 65536, maximum_solves = 65536;
  // Result header + rows; excludes borrowed inputs, stack, allocator and RSS.
  std::size_t maximum_native_bytes = 1ull << 30;
};
struct HydrogenBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::size_t solves = 0;
  std::vector<HydrogenRow> rows;
};
// Charge neutrality, ideal nonrelativistic Maxwell-Boltzmann species, ground
// state only, blackbody zero photon chemical potential; fixed translational
// mass ratio m_p/m_H=1. Spin counts e/p/H=2/1/2 give prefactor 1.
// Finite 1<=T<=1e5 K and n_H>0, with conservative n_H/n_Qe<=1e-3 admission;
// n_Qe=(2*pi*m_e*kB*T/h^2)^1.5. This does not certify LTE or excited-state
// approximation accuracy. No cosmic history, helium or drag prediction.
// Ordered coarse batch, exact source scalar bits, independent requested groups.
// Positive fractions must meet 1e-13 relative arithmetic/cast admission; never
// replace unrepresentable fractions by zero. A complementary rounded value may
// equal one while its independently stored tiny positive partner remains >0.
// Native allocation failure propagates ordinary std::bad_alloc, with RAII.
HydrogenBatch evaluate_hydrogen_equilibrium(std::span<const HydrogenState>,
                                            HydrogenPolicy = {});
} // namespace irred::atomic
