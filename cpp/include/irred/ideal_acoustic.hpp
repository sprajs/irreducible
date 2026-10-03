#pragma once
#include "irred/thermal_neutrino.hpp"
#include <array>
#include <string>

namespace irred::cosmology {
// Array/bit order is stable. Every coordinate is dimensionless per asymptotic
// unit zeta; theta/k is the signed MB Doppler coefficient, not theta itself.
inline constexpr unsigned acoustic_comoving_cdm = 1;
inline constexpr unsigned acoustic_comoving_baryon = 2;
inline constexpr unsigned acoustic_intrinsic_temperature = 4;
inline constexpr unsigned acoustic_metric_phi = 8;
inline constexpr unsigned acoustic_newtonian_cdm = 16;
inline constexpr unsigned acoustic_newtonian_baryon = 32;
inline constexpr unsigned acoustic_doppler_baryon = 64;
inline constexpr unsigned acoustic_doppler_cdm = 128;
inline constexpr unsigned acoustic_metric_log_derivative = 256;
inline constexpr unsigned acoustic_all_outputs = 511;
inline constexpr std::size_t acoustic_output_count = 9;
struct IdealAcousticRequest {
  ThermalPhysicalModel model;
  double initial_scale_factor = 1e-10;
  std::string source_origin;
};
struct IdealAcousticPolicy {
  double absolute_tolerance = 1e-8, relative_tolerance = 3e-5;
  double maximum_log_step = .02, maximum_constraint_residual = 1e-6;
  std::size_t maximum_wavenumbers = 2, maximum_samples = 4096;
  // Original defaults remain. The augmented caller may explicitly request at
  // most 2M per k/4M per batch for nominal plus six response-vector owners.
  std::size_t maximum_rhs_per_wavenumber = 1000000;
  std::size_t maximum_rhs_batch = 2000000;
  std::size_t maximum_background_evaluations = 3000000;
  std::size_t maximum_age_evaluations = 1000000;
  std::size_t maximum_state_element_writes = 30000000;
  std::size_t maximum_diagnostic_evaluations = 3000000;
  std::size_t maximum_native_bytes = 16 * 1024 * 1024;
};
struct IdealAcousticWork {
  std::size_t physical_mappings = 0, background_preparations = 0;
  std::size_t attempts = 0, rhs_evaluations = 0;
  std::size_t background_evaluations = 0, age_evaluations = 0;
  std::size_t state_element_writes = 0, endpoint_evaluations = 0;
  std::size_t diagnostic_evaluations = 0, output_evaluations = 0;
};
struct IdealAcousticError {
  double time_refinement = 0, initial_refinement = 0;
  std::optional<double> common_source_background_age;
  std::optional<double> arithmetic_storage_constraint;
  std::optional<double> absolute_error_estimate;
};
struct IdealAcousticValue {
  numerics::Status status = numerics::Status::invalid_input;
  // Computed is the actual finite numerical witness, even when its source/error
  // gate refuses. It is never an admitted prediction without status/value.
  std::optional<long double> computed;
  std::optional<double> value;
  IdealAcousticError error;
};
// Operation-owned first-refusal IDs. Stage0 means no RK stage; stages1..4
// identify the attempted RK input. These are diagnostics, not output IDs.
enum class IdealAcousticFailureStage : unsigned {
  initial_age = 1, initial_background, initial_loading, initial_domain,
  initial_state_storage, initial_projection, initial_clock, initial_frame,
  initial_source_bounds, initial_response, sample_storage, step_clock,
  step_background, rk_storage, rk_nominal_rhs, rk_frame, rk_response_rhs,
  rk_combine, endpoint_combine, committed_state, hamiltonian_constraint,
  phase_response, phase_domain, endpoint_storage, radius_diagonal,
  radius_forcing, radius_predictor, radius_local_shadow, radius_corrector,
  radius_endpoint_bridge, radius_endpoint_sum
};
// These identify the six stored evolution coordinates, not requested output
// IDs. Channel values retain the literal response-vector offsets.
enum class IdealAcousticStateCoordinate : unsigned {
  comoving_cdm = 0, cdm_photon_entropy, scaled_common_minus_cdm_velocity,
  scaled_cdm_velocity, metric_phi, conformal_age
};
enum class IdealAcousticRadiusChannel : unsigned {
  source = 24, arithmetic = 30
};
struct IdealAcousticRadiusFailure {
  IdealAcousticStateCoordinate coordinate = IdealAcousticStateCoordinate::comoving_cdm;
  IdealAcousticRadiusChannel channel = IdealAcousticRadiusChannel::source;
  // Borrowed scalars from the failing guard, without a repeat calculation.
  // Provisional includes the comparison-sum rounding allowance. The local
  // assembly allowance exists only for ARITHMETIC. Nonfinite scalars are absent;
  // finite negative/subnormal scalars remain visible for diagnosis.
  std::optional<long double> provisional_radius;
  std::optional<long double> local_assembly_allowance;
  std::optional<long double> completed_radius;
  numerics::Status arithmetic_status = numerics::Status::invalid_input;
  bool completed_finite = false, completed_normal_or_zero = false;
  bool completed_nonnegative = false;
};
struct IdealAcousticFailure {
  IdealAcousticFailureStage stage = IdealAcousticFailureStage::initial_age;
  unsigned rk_stage = 0;
  // Last fully committed nominal state, before a failing partial assembly.
  // Unavailable before the initial projected nominal state exists.
  std::optional<long double> last_committed_scale_factor;
  // Existing coefficient/input epoch used by the failing operation; no query
  // is made to fill this witness. Clock/storage-only failures may lack it.
  std::optional<long double> attempted_coefficient_scale_factor;
  // Present only at the existing completed-radius guard in an RK combination.
  // Original attempt status and first-refusal stage remain authoritative.
  std::optional<IdealAcousticRadiusFailure> radius_assembly;
};
struct IdealAcousticAttempt {
  numerics::Status status = numerics::Status::invalid_input;
  double initial_scale_factor = 0, maximum_log_step = 0;
  std::optional<std::array<long double, 5>> unprojected_initial_state;
  std::optional<std::array<long double, 5>> projected_initial_state;
  std::optional<long double> initial_delta_projection;
  // Successful attempts have no failure. Original status remains authoritative.
  std::optional<IdealAcousticFailure> failure;
  long double maximum_normalized_hamiltonian_residual = 0;
  long double maximum_absolute_hamiltonian_residual = 0;
  long double maximum_relative_reduced_constraint = 0;
  long double maximum_direct_hamiltonian_residual = 0;
  long double maximum_hamiltonian_assembly_discrepancy = 0;
  long double maximum_absolute_closure_defect = 0;
  long double maximum_absolute_acceleration_defect = 0;
  IdealAcousticWork work;
};
struct IdealAcousticTrajectory {
  std::size_t original_k_index = 0;
  double wavenumber_mpc_inverse = 0;
  std::array<IdealAcousticAttempt, 5> attempts;
  std::size_t attempts_started = 0;
  IdealAcousticWork work;
};
struct IdealAcousticEpoch {
  numerics::Status status = numerics::Status::invalid_input;
  long double state_scale_factor = 0, hcal_mpc_inverse = 0;
  long double conformal_age_mpc = 0, endpoint_scale_factor_error = 0;
  // Missing source/derivative ownership stays unavailable, never zero.
  std::optional<long double> age_error_mpc, derivative_consistency_estimate;
};
struct IdealAcousticRow {
  std::size_t original_k_index = 0, original_a_index = 0;
  double wavenumber_mpc_inverse = 0, requested_scale_factor = 0;
  std::array<IdealAcousticValue, acoustic_output_count> outputs;
  IdealAcousticEpoch epoch;
};
struct IdealAcousticBatch {
  numerics::Status status = numerics::Status::invalid_input;
  numerics::Status shared_dependency_status = numerics::Status::invalid_input;
  unsigned requested_outputs = 0;
  std::vector<IdealAcousticRow> rows; // k-major, a-minor; no sorting/deduplication.
  std::vector<IdealAcousticTrajectory> trajectories;
  IdealAcousticWork preparation_work, evaluation_work;
  std::optional<std::size_t> peak_owned_payload_bound;
};
class IdealAcousticTransfer {
public:
  IdealAcousticTransfer() = default;
  IdealAcousticTransfer(const IdealAcousticTransfer &) = delete;
  IdealAcousticTransfer &operator=(const IdealAcousticTransfer &) = delete;
  IdealAcousticTransfer(IdealAcousticTransfer &&) noexcept;
  IdealAcousticTransfer &operator=(IdealAcousticTransfer &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const IdealAcousticRequest *source() const noexcept { return source_ ? &*source_ : nullptr; }
  const ThermalPhysicalMapping *physical_mapping() const noexcept { return mapping_ ? &*mapping_ : nullptr; }
  const ThermalBackground *background() const noexcept { return background_ ? &*background_ : nullptr; }
  IdealAcousticWork preparation_work() const noexcept { return preparation_work_; }
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  // a must be strictly increasing and lie in [source ai,.01]. All five solves
  // retain the same requested endpoints; original k order and repeats survive.
  IdealAcousticBatch evaluate(std::span<const double> k_mpc_inverse,
                             std::span<const double> scale_factors,
                             unsigned outputs, IdealAcousticPolicy = {}) const;
private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<IdealAcousticRequest> source_;
  std::optional<ThermalPhysicalMapping> mapping_;
  std::optional<ThermalBackground> background_;
  long double loading_ratio_ = 0, loading_error_ = 0;
  long double loading_numerator_ = 0, loading_denominator_ = 0;
  long double retained_lambda_ = 0, lambda_getter_signed_loss_ = 0;
  // Facts from the single preparation capture, retained for the pending common
  // diagnostic with the original mapper witnesses. They avoid a new capture,
  // remap or reporting getter call to recreate that diagnostic's ancestry.
  long double retained_critical_density_ = 0, retained_species_today_ = 0;
  long double retained_species_normalization_error_ = 0;
  double retained_lambda_emitted_ = 0;
  ThermalMomentumMethod retained_momentum_method_ = ThermalMomentumMethod::direct_adaptive;
  IdealAcousticWork preparation_work_;
  friend IdealAcousticTransfer prepare_ideal_acoustic(IdealAcousticRequest &&, IdealAcousticPolicy);
};
IdealAcousticTransfer prepare_ideal_acoustic(IdealAcousticRequest &&,
                                           IdealAcousticPolicy = {});
inline constexpr std::string_view ideal_acoustic_model_id =
    "GR/flat-ideal-infinite-coupling-photon-cold-baryon-CDM-smooth-Lambda-unit-zeta/v1";
inline constexpr std::string_view ideal_acoustic_method_id =
    "central-RK4/radius-PC2/coupled-local-RK4-impulses/v2";
inline constexpr std::string_view ideal_acoustic_arithmetic_id =
    "strict-wide-conditional-64eps-clock-common-thermal-state-nine-signed-readouts/v2";
} // namespace irred::cosmology
