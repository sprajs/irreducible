#pragma once
#include "irred/hydrogen_helium_equilibrium.hpp"
#include <string>
namespace irred::cosmology {
inline constexpr std::string_view baryon_abundance_model_id =
    "supplied-H1-He4-neutral-effective-mass-abundance/v1";
inline constexpr std::string_view baryon_abundance_method_id =
    "shared-H100-critical-mass-density-ordered-a-cubed/v1";
inline constexpr std::string_view baryon_abundance_assets_id =
    "SI2019-IAU2015-Mpc-CODATA2018-G-caller-neutral-masses";
inline constexpr std::string_view baryon_abundance_lte_model_id =
    "supplied-H1-He4-abundance-shared-electron-ground-state-LTE/v1";
struct BaryonAbundanceSource {
  // omega_b=Omega_b*h^2, h=H0/100; no H0 argument is needed here.
  double physical_baryon_density = 0, helium4_mass_fraction = 0;
  // Supplied positive normal neutral ground-state effective masses, kg.
  double hydrogen1_effective_mass_kg = 0, helium4_effective_mass_kg = 0;
  std::string source_origin, mass_origin;
};
struct BaryonAbundancePolicy {
  std::size_t maximum_rows = 65536, maximum_native_bytes = 1ull << 30;
};
struct BaryonAbundanceValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double relative_arithmetic_estimate = 0;
};
struct BaryonAbundanceRow {
  double scale_factor = 0;
  numerics::Status admission_status = numerics::Status::invalid_input;
  // Physical total nuclei per m^3, counting all ionization stages.
  BaryonAbundanceValue hydrogen_nuclei, helium_nuclei;
};
struct BaryonAbundanceBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::vector<BaryonAbundanceRow> rows;
};
struct BaryonAbundanceLteQuery {
  double scale_factor = 0, temperature_kelvin = 0;
};
struct BaryonAbundanceLteRow {
  BaryonAbundanceLteQuery source;
  BaryonAbundanceRow nuclei;
  atomic::HydrogenHeliumRow equilibrium;
};
struct BaryonAbundanceLteBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::string temperature_origin;
  std::size_t solves = 0, root_iterations = 0, charge_evaluations = 0;
  std::vector<BaryonAbundanceLteRow> rows;
};
class BaryonAbundance {
public:
  BaryonAbundance() = default;
  BaryonAbundance(const BaryonAbundance &) = default;
  BaryonAbundance &operator=(const BaryonAbundance &);
  BaryonAbundance(BaryonAbundance &&) noexcept;
  BaryonAbundance &operator=(BaryonAbundance &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const std::optional<BaryonAbundanceSource> &source() const noexcept { return source_; }
  BaryonAbundanceBatch evaluate(std::span<const double>, BaryonAbundancePolicy = {}) const;
  BaryonAbundanceLteBatch evaluate_equilibrium(
      std::span<const BaryonAbundanceLteQuery>, std::string_view temperature_origin,
      atomic::HydrogenHeliumPolicy = {}, BaryonAbundancePolicy = {}) const;
private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<BaryonAbundanceSource> source_;
  long double rho0_ = 0, rho_error_ = 0;
  BaryonAbundanceRow map(double) const;
  std::optional<std::size_t> payload() const noexcept;
  friend BaryonAbundance prepare_baryon_abundance(const BaryonAbundanceSource &, BaryonAbundancePolicy);
};
// Supplied Y is the He4 share of the declared fixed neutral-effective mass
// density, not a number fraction. rho0=rho_crit(H100)*omega_b;
// nH=rho0*(1-Y)/(a^3*mH), nHe=rho0*Y/(a^3*mHe).
// No BBN prediction, evolving binding mass, kinetics or expansion feedback.
// Mapping supports exact zero endpoints; LTE requires both species present.
// Strict nearest wide arithmetic (>=64 bits), mapping <=1e-15 relative
// diagnostic; combined LTE <=1e-13. Physical uncertainty is excluded.
// Explicit payload bounds exclude allocator bookkeeping/RSS and borrowed inputs.
// Allocation failure propagates bad_alloc with RAII cleanup.
BaryonAbundance prepare_baryon_abundance(const BaryonAbundanceSource &, BaryonAbundancePolicy = {});
} // namespace irred::cosmology
