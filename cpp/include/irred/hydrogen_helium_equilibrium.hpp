#pragma once
#include "irred/numerics.hpp"
#include <string_view>
namespace irred::atomic {
// Fixed ASD5.12 central scalar facts for 4He. Physical asset uncertainty is
// excluded from numerical diagnostics; HeII's ASD energy is theoretical.
inline constexpr long double helium_first_ionization_energy_ev = 24.587389011L;
inline constexpr long double helium_first_ionization_uncertainty_ev = 2.5e-8L;
inline constexpr long double helium_second_ionization_energy_ev =
    54.4177655282L;
inline constexpr long double helium_second_ionization_uncertainty_ev = 1e-9L;
inline constexpr std::string_view hydrogen_helium_equilibrium_model_id =
    "homogeneous_ground_state_hydrogen_helium_saha_supplied_temperature_nuclei_"
    "densities";
inline constexpr std::string_view hydrogen_helium_equilibrium_assets_id =
    "SI2019-CODATA2022-me-NIST-ASD5.12-H1-He4-fixed";
inline constexpr std::string_view hydrogen_helium_equilibrium_method_id =
    "log_neutrality_safeguarded_newton_separate_positive_stage_fractions";
inline constexpr std::string_view hydrogen_helium_equilibrium_arithmetic_id =
    "longdouble64-cpu-nearest-strict-fp-v1";
enum HydrogenHeliumOutput : unsigned {
  hydrogen_neutral_fraction = 1,
  hydrogen_ionized_fraction = 2,
  helium_neutral_fraction = 4,
  helium_singly_ionized_fraction = 8,
  helium_doubly_ionized_fraction = 16,
  electron_number_density = 32
};
struct HydrogenHeliumState {
  double temperature_kelvin = 0;
  // Physical total nuclei densities (all stages), independently supplied m^-3.
  double hydrogen_nuclei_per_cubic_metre = 0;
  double helium_nuclei_per_cubic_metre = 0;
};
struct HydrogenHeliumValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  // Empirical arithmetic/libm/root/cast diagnostic, not atomic/model error.
  double relative_arithmetic_estimate = 0;
};
struct HydrogenHeliumRow {
  HydrogenHeliumState source;
  numerics::Status admission_status = numerics::Status::invalid_input;
  std::size_t solves = 0, root_iterations = 0, charge_evaluations = 0;
  HydrogenHeliumValue hydrogen_neutral, hydrogen_ionized, helium_neutral,
      helium_singly_ionized, helium_doubly_ionized;
  // Physical free electrons m^-3, charge neutral with both supplied species.
  HydrogenHeliumValue electron_density;
};
struct HydrogenHeliumPolicy {
  unsigned requested_outputs = 63;
  std::size_t maximum_rows = 65536, maximum_solves = 65536;
  std::size_t maximum_root_iterations = 256;
  std::size_t maximum_charge_evaluations = 4000000;
  // Result header + rows; excludes borrowed inputs, stack, allocator and RSS.
  std::size_t maximum_native_bytes = 1ull << 30;
};
struct HydrogenHeliumBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::size_t solves = 0, root_iterations = 0, charge_evaluations = 0;
  std::vector<HydrogenHeliumRow> rows;
};
// Homogeneous ground-state H/4He LTE with zero photon chemical potential,
// dilute nonrelativistic Maxwell-Boltzmann species and adjacent translational
// mass ratios approximated by one. Qe is the unchanged hydrogen-owned equation;
// stage Saha factors H/HeI/HeII are Qe/4Qe/Qe times their Boltzmann factors.
// Shared electron neutrality; no abundance, cosmic density or kinetic mapping.
// Finite 1<=T<=1e5 K, nH,nHe>0; conservative (nH+2nHe)/Qe<=1e-3.
// Six independently requested positive outputs, <=1e-13 relative diagnostic.
// Failed/unrequested partners are never zeros; dominant fractions may round1.
// Source scalar bits/order owned by rows. Bounds precede scans/allocation.
// Allocation failure propagates std::bad_alloc with ordinary RAII cleanup.
HydrogenHeliumBatch
    evaluate_hydrogen_helium_equilibrium(std::span<const HydrogenHeliumState>,
                                         HydrogenHeliumPolicy = {});
} // namespace irred::atomic
