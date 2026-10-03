#pragma once
#include "irred/gaussian_box.hpp"
#include <memory>

namespace irred::statistics {
struct BoxHeldoutMetadata {
  std::string source_contract_id, candidate_identity, conditioning_identity;
  std::string original_row_lineage, event_lineage, calibration_dependence_identity;
  std::string heldout_unit, heldout_covariance_unit, heldout_measure;
  // Explicit original47->active46 mapping and literal point-mass witness.
  std::vector<std::string> ordered_original_parameter_ids;
  std::vector<std::size_t> active_original_parameter_indices;
  std::size_t fixed_original_parameter_index = 0;
  double fixed_original_parameter_value = 0;
  std::string fixed_coordinate_provenance;
};
struct BoxHeldoutPolicy {
  BoxPolicy box;
  std::size_t maximum_training_vectors = 1, maximum_candidate_values = 1,
      maximum_requests = 1;
  // Required caller-supplied frozen caps; zero refuses. Upper charges, not
  // measured instruction counts; overflow and failed attempts are retained.
  std::size_t maximum_preparation_work_units = 0, maximum_evaluation_work_units = 0;
  double maximum_conditional_log_density_width = 3e-9;
};
enum class BoxHeldoutStage {
  unassessed, source_validation, block_extraction, training_factor,
  cross_whitening, schur_admission, training_qr, update_qr,
  training_normalization, joint_normalization, ratio_composition, complete
};
// Logical loop-body/element-copy/helper-entry upper charges, not FLOPs or CPU
// instructions. Every attempted primitive is precharged before its first loop.
struct BoxHeldoutWork {
  std::size_t factor_attempts = 0, factors_completed = 0,
      whitening_attempts = 0, whitenings_completed = 0,
      training_qr_attempts = 0, training_qr_completed = 0,
      update_qr_attempts = 0, update_qr_completed = 0,
      training_attempts = 0, training_completed = 0,
      joint_attempts = 0, joint_completed = 0,
      variance_screen_attempts = 0, variance_screens_completed = 0,
      tail_bound_attempts = 0, tail_bounds_completed = 0,
      normalization_attempts = 0, normalizations_completed = 0,
      composition_attempts = 0, compositions_completed = 0,
      charged_work_units = 0;
};
struct BoxHeldoutPreparation {
  DensityStatus status = DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  BoxHeldoutStage stage = BoxHeldoutStage::unassessed;
  bool cross_whitening_available = false, schur_available = false,
      schur_sqrt_available = false, appended_response_available = false,
      training_qr_available = false, update_qr_available = false,
      covariance_determinants_available = false, retained_frames_available = false;
  // Encloses sqrt of the separately reported binary64 schur_variance;
  // internal whitening uses the wide Schur proposal, screened separately.
  BoxInterval schur_sqrt_enclosure;
  double schur_variance = 0, schur_cancellation_ratio = 0,
      schur_sensitivity_estimate = 0, cross_whitening_backward_residual = 0,
      cross_whitening_rounding_estimate = 0,
      appended_projection_rounding_estimate = 0,
      training_preparation_sensitivity = 0, training_triangular_condition = 0, training_transpose_condition = 0,
      update_triangular_condition = 0, update_transpose_condition = 0,
      full_source_covariance_log_determinant = 0,
      training_covariance_log_determinant = 0;
  std::optional<std::size_t> peak_payload_bound, retained_payload_bound;
  // Present only if this failure actually produced a scalar witness. May be
  // nonfinite; external transport must tag it rather than invent a finite zero.
  std::optional<long double> refusal_witness;
  std::optional<std::size_t> refusal_index;
  BoxHeldoutWork work;
};
struct BoxHeldoutRequest {
  std::size_t training_vector_index = 0, candidate_value_index = 0;
};
struct BoxHeldoutDensity {
  DensityStatus status = DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  BoxHeldoutStage stage = BoxHeldoutStage::unassessed;
  bool density_available = false, scalar_conditional_constant_available = false,
      cross_projection_diagnostics_available = false,
      offset_diagnostics_available = false, update_diagnostics_available = false;
  BoxInterval log_density, scalar_conditional_log_normalization;
  // Shared normalization helper fills explicit availability; quantile/CDF
  // fields are unavailable and their work is zero for this consumer.
  GaussianBoxResult joint_normalization;
  double cross_projection_rounding_estimate = 0,
      offset_subtraction_rounding_estimate = 0,
      update_rhs_rounding_estimate = 0;
  std::optional<long double> refusal_witness;
  std::optional<std::size_t> refusal_index;
  BoxHeldoutWork work;
};
struct BoxHeldoutBatch {
  // finite only admits the batch shape; each training/result status still owns
  // its availability. Shared admission refusals publish no usable scalar.
  DensityStatus status = DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  // False on a partial outer allocation prefix. Individual arrays may then
  // have different sizes; only present availability flags qualify scalars.
  bool output_layout_available = false;
  std::vector<GaussianBoxResult> training_normalizations;
  std::vector<double> training_offset_subtraction_rounding_estimates;
  // Absolute cached dot-product error estimate, not a certificate. Kept
  // outside the exact 2p+2 wide cache and used again by each queried subtraction.
  std::vector<double> training_cross_projection_error_estimates;
  std::vector<std::uint8_t> training_cross_projection_diagnostics_available;
  std::vector<std::optional<long double>> training_refusal_witnesses;
  std::vector<std::uint8_t> training_offset_diagnostics_available;
  std::vector<BoxHeldoutDensity> densities;
  std::vector<BoxHeldoutRequest> requests;
  BoxHeldoutWork work;
  std::optional<std::size_t> scratch_output_payload_bound;
  const char *method_id = "retained-training-qr-rankone-finitebox-evidence-ratio/v1";
  const char *law_id = "finite-uniform-box-fixed-design-correlated-one-row-predictive/v1";
  const char *enclosure_scope = "conditional-on-two-reported-completions;not-original-input-certificate/v1";
};

class GaussianBoxHeldout {
public:
  GaussianBoxHeldout() noexcept;
  ~GaussianBoxHeldout(); // out-of-line because Impl is incomplete here
  GaussianBoxHeldout(const GaussianBoxHeldout &) = delete;
  GaussianBoxHeldout &operator=(const GaussianBoxHeldout &) = delete;
  GaussianBoxHeldout(GaussianBoxHeldout &&) noexcept;
  GaussianBoxHeldout &operator=(GaussianBoxHeldout &&) noexcept;
  // X contains all original rows and the active columns in the declared
  // original-axis map. The separately witnessed fixed coordinate is literal0;
  // the untouched original full design remains external input lineage.
  // Failure leaves full_source usable. Success consumes it only after every
  // block/factor/update/admission/budget gate. Arrays are borrowed at setup.
  static GaussianBoxHeldout prepare(Gaussian &&full_source,
      std::span<const double> full_original_row_major_design,
      std::span<const double> full_original_offsets,
      std::span<const std::string> full_original_row_ids,
      DesignMetadata, BoxSupport, std::size_t heldout_original_row_index,
      BoxHeldoutMetadata, BoxHeldoutPolicy);
  static std::optional<std::size_t> preparation_payload_bound(
      const Gaussian &, const DesignMetadata &, const BoxSupport &,
      const BoxHeldoutMetadata &) noexcept;
  // These total upper envelopes include all source/factor/QR/normalizer paths,
  // including failed prefixes. Callers freeze an explicit cap >= the returned
  // total; zero caps refuse. They are not runtime instructions or new budgets.
  static std::optional<std::size_t> preparation_work_bound(
      const Gaussian &, const DesignMetadata &, const BoxSupport &,
      const BoxHeldoutMetadata &) noexcept;
  static std::optional<std::size_t> evaluation_work_bound(
      std::size_t full_dimension,std::size_t parameter_count,
      std::size_t training_vectors,std::size_t candidate_values,std::size_t requests,
      std::size_t training_identity_bytes) noexcept;
  std::optional<std::size_t> evaluation_work_bound(
      std::size_t training_vectors, std::size_t candidate_values,std::size_t requests) const noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  std::optional<std::size_t> evaluation_payload_bound(
      std::size_t training_vectors, std::size_t requests) const noexcept;
  DensityStatus status() const noexcept;
  const BoxHeldoutPreparation &preparation() const noexcept;
  std::span<const std::string> original_row_ids() const noexcept;
  std::span<const std::string> training_row_ids() const noexcept;
  const BoxHeldoutMetadata &metadata() const noexcept;
  // Each training vector has exactly n-1 values in training_row_ids order.
  // Pool/candidate/request order is retained. The withheld value is never a
  // member of a training pool. Repeated requests do not refactor C or X.
  BoxHeldoutBatch evaluate(std::span<const double> pooled_training_values,
      std::span<const std::string> training_row_ids,
      std::span<const double> pooled_candidate_values,
      std::span<const BoxHeldoutRequest>, BoxHeldoutPolicy) const;
private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  // Explicit moves reset source receipt/accessors, including failed prepare.
  BoxHeldoutPreparation preparation_;
};
} // namespace irred::statistics
