#pragma once
#include "irred/baryon_abundance.hpp"
#include <memory>
namespace irred::cosmology {
inline constexpr std::string_view baryon_abundance_density_law_model_id =
    "finite-joint-supplied-H1-He4-neutral-mass-density-law/v1";
inline constexpr std::string_view baryon_abundance_lte_law_model_id =
    "finite-joint-supplied-H1-He4-abundance-mass-matter-temperature-LTE-law/v1";
inline constexpr std::string_view baryon_abundance_law_method_id =
    "retained-abundance-owners-wide-normalization-centered-population-moments/v1";
inline constexpr std::string_view baryon_abundance_law_arithmetic_id =
    "longdouble64-cpu-nearest-strict-fp-v1";
struct BaryonAbundanceLawRow {
  std::string_view id;
  double scale_factor = 0;
};
struct BaryonAbundanceLawState {
  std::string_view id;
  double relative_mass = 0;
  BaryonAbundanceSource abundance;
  // Entire joint state's supplied matter temperatures, in declared row order.
  std::span<const double> temperature_kelvin;
  std::string_view temperature_origin;
};
struct BaryonAbundanceLawInput {
  std::span<const BaryonAbundanceLawRow> rows;
  std::span<const BaryonAbundanceLawState> states;
  std::string_view distribution_origin, dependence_origin;
};
struct BaryonAbundanceLawPolicy {
  std::size_t maximum_states = 1024, maximum_rows = 64;
  std::size_t maximum_state_rows = 65536, maximum_axes = 384;
  std::size_t maximum_covariance_products = 1ull << 27;
  // Solves and charge evaluations are cumulative across all child calls,
  // including refused solves; maximum_root_iterations is per solve.
  std::size_t maximum_solves = 65536, maximum_root_iterations = 256;
  std::size_t maximum_charge_evaluations = 4000000;
  std::size_t maximum_native_bytes = 1ull << 30;
  // Empirical absolute diagnostics compared to |mean| and sqrt(Cii*Cjj).
  // The named contract may be tightened, not relaxed by these fields.
  double mean_relative_sensitivity = 1e-12;
  double covariance_relative_sensitivity = 1e-8;
};
struct BaryonAbundanceLawOwnedRow {
  std::string id;
  double scale_factor = 0;
};
struct BaryonAbundanceLawRetainedState {
  std::string id;
  double relative_mass = 0;
  // Numerical coordinates: no binary64 probability bottleneck or lost state.
  long double probability = 0, probability_relative_arithmetic_estimate = 0;
  BaryonAbundance abundance;
  std::vector<double> temperature_kelvin;
  std::string temperature_origin;
};
struct BaryonAbundanceLawSource {
  std::string distribution_origin, dependence_origin;
  std::vector<BaryonAbundanceLawOwnedRow> rows;
  std::vector<BaryonAbundanceLawRetainedState> states;
  std::size_t state_preparations = 0, normalization_terms = 0;
  // Conservative copied-string envelope checked before preparation allocation.
  std::size_t preparation_native_payload_bound = 0;
  // Explicit retained payload; excludes allocator bookkeeping/control block/RSS.
  std::size_t retained_native_payload_bytes = 0;
};
enum class BaryonAbundanceLawOutput {
  hydrogen_nuclei, helium_nuclei,
  hydrogen_neutral, hydrogen_ionized, helium_neutral,
  helium_singly_ionized, helium_doubly_ionized, electron_density
};
struct BaryonAbundanceLawAxis {
  std::size_t row_index = 0;
  BaryonAbundanceLawOutput output = BaryonAbundanceLawOutput::hydrogen_nuclei;
};
struct BaryonAbundanceLawMoments {
  std::vector<double> mean, covariance; // Full row-major population covariance.
  std::vector<double> mean_empirical_sensitivity,
      covariance_empirical_sensitivity; // Absolute, corresponding product units.
};
struct BaryonAbundanceLawWork {
  std::size_t mapping_rows = 0, solves = 0, root_iterations = 0,
      charge_evaluations = 0, covariance_products = 0;
  std::size_t combined_native_payload_bound = 0;
};
struct BaryonAbundanceDensityLawResult {
  numerics::Status status = numerics::Status::invalid_input;
  std::shared_ptr<const BaryonAbundanceLawSource> source;
  std::vector<BaryonAbundanceLawAxis> axes; // Row, then H before He.
  std::vector<BaryonAbundanceBatch> attempts; // State order; nothing dropped.
  BaryonAbundanceLawWork work;
  std::optional<BaryonAbundanceLawMoments> moments;
};
struct BaryonAbundanceLteLawResult {
  numerics::Status status = numerics::Status::invalid_input;
  std::shared_ptr<const BaryonAbundanceLawSource> source;
  std::vector<BaryonAbundanceLawAxis> axes; // Row, then ascending native mask bit.
  std::vector<BaryonAbundanceLteBatch> attempts; // State order, including failures.
  BaryonAbundanceLawWork work;
  std::optional<BaryonAbundanceLawMoments> moments;
};
class BaryonAbundanceLaw {
public:
  BaryonAbundanceLaw() = default;
  BaryonAbundanceLaw(const BaryonAbundanceLaw &) = default;
  BaryonAbundanceLaw &operator=(const BaryonAbundanceLaw &) = default;
  BaryonAbundanceLaw(BaryonAbundanceLaw &&) noexcept;
  BaryonAbundanceLaw &operator=(BaryonAbundanceLaw &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const std::shared_ptr<const BaryonAbundanceLawSource> &source() const noexcept {
    return source_;
  }
  BaryonAbundanceDensityLawResult evaluate_density(BaryonAbundanceLawPolicy = {}) const;
  BaryonAbundanceLteLawResult evaluate_equilibrium(
      unsigned requested_outputs = 63, BaryonAbundanceLawPolicy = {}) const;
private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::shared_ptr<const BaryonAbundanceLawSource> source_;
  friend BaryonAbundanceLaw prepare_baryon_abundance_law(
      const BaryonAbundanceLawInput &, BaryonAbundanceLawPolicy);
};
// Finite FULL JOINT supplied omega_b/Y/neutral-effective-mass/matter-T law.
// Positive normal binary64 masses and T; fixed 0<a<=1; all source roles/IDs owned.
// Normalize once wide. Density endpoints remain valid; LTE keeps its own domain.
// Existing abundance/LTE equations are the sole physical owners. G and atomic
// assets remain fixed; no BBN, kinetics, posterior or physical qualification.
// Required state/output failures preserve attempts but withhold all moments.
// Mathematical constant-input witnesses alone admit zero variance; equal rounded
// outputs of varying inputs are not witnesses. No jitter or renormalized subset.
// Density/electrons: m^-3; fractions: dimensionless; covariance: product units.
// Combined explicit payload includes retained law, outputs, scratch and peak
// child temporaries; excludes borrowed inputs, stack, allocator/control block/RSS.
// std::bad_alloc propagates with ordinary RAII cleanup.
BaryonAbundanceLaw prepare_baryon_abundance_law(
    const BaryonAbundanceLawInput &, BaryonAbundanceLawPolicy = {});
} // namespace irred::cosmology
