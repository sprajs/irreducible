#pragma once
#include "irred/gaussian_posterior.hpp"
namespace irred::statistics {
struct PredictiveMetadata {
  std::vector<std::string> ordered_parameter_ids, parameter_units,
      response_units, training_event_ids, future_event_ids;
  std::string future_unit, future_covariance_unit, future_measure,
      response_identity, conditioning_identity, conditional_noise_identity,
      dependence_identity;
  bool future_noise_independence_declared = false;
  bool noise_conditional_on_parameters_declared = false;
};
struct PredictivePolicy {
  std::size_t maximum_elements = 1000000;
  std::size_t maximum_payload_bytes = 256 * 1024 * 1024;
  std::size_t maximum_work_units = 100000000;
  double maximum_forward_sensitivity = 1e-10;
};
// Fixed training r and response A: y*|r ~ N(A mu,R*+A V A^T).
// The explicit independent future noise is conditional on the same beta.
// This is a derived predictive law, not an additional observed measurement.
class GaussianPredictive {
public:
  GaussianPredictive() = default;
  GaussianPredictive(const GaussianPredictive &) = delete;
  GaussianPredictive &operator=(const GaussianPredictive &) = delete;
  GaussianPredictive(GaussianPredictive &&) noexcept;
  GaussianPredictive &operator=(GaussianPredictive &&) noexcept;
  // Both immutable inputs are borrowed only during preparation and remain
  // usable on success/failure. No factor is copied from either input.
  static GaussianPredictive
  prepare(const GaussianPosterior &, std::span<const double> training_values,
          std::span<const std::string> training_row_ids,
          const Gaussian &future_noise,
          std::span<const double> row_major_response,
          const PredictiveMetadata &, PredictivePolicy = {});
  // Includes both already-retained borrowed owners once, all simultaneous
  // candidate storage and metadata. Excludes borrowed arrays/allocator/RSS.
  static std::optional<std::size_t>
  preparation_payload_bound(const GaussianPosterior &, const Gaussian &,
                            const PredictiveMetadata &) noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  // Scratch only, excluding the retained owner and borrowed future vector.
  std::optional<std::size_t> evaluation_payload_bound() const noexcept;
  DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  numerics::Arithmetic arithmetic() const noexcept {
    return predictive_ ? predictive_->arithmetic()
                       : numerics::Arithmetic::binary64_legacy_v1;
  }
  const PredictiveMetadata &metadata() const noexcept { return metadata_; }
  const Metadata &training_source_metadata() const noexcept {
    return training_source_;
  }
  const Metadata &future_noise_metadata() const noexcept {
    return future_noise_;
  }
  const std::string &prior_identity() const noexcept { return prior_identity_; }
  const std::string &training_design_identity() const noexcept {
    return training_design_identity_;
  }
  const std::string &parameter_measure() const noexcept {
    return parameter_measure_;
  }
  const std::string &prior_dependence_identity() const noexcept {
    return prior_dependence_identity_;
  }
  std::span<const std::string> shared_nuisance_ids() const noexcept {
    return shared_nuisance_ids_;
  }
  std::span<const double> training_values() const noexcept {
    return training_values_;
  }
  std::span<const double> response() const noexcept { return response_; }
  // These views borrow this owner; moves/destruction invalidate them.
  std::span<const double> mean() const noexcept { return mean_; }
  std::span<const double> mean_absolute_error_estimates() const noexcept {
    return mean_error_;
  }
  std::span<const double> covariance() const noexcept {
    return predictive_ ? predictive_->covariance() : std::span<const double>{};
  }
  std::span<const double> covariance_absolute_error_estimates() const noexcept {
    return covariance_error_;
  }
  double covariance_relative_error_estimate() const noexcept {
    return relative_covariance_error_;
  }
  std::size_t preparation_work_units() const noexcept {
    return preparation_work_units_;
  }
  GaussianResult log_density(std::span<const double> future_values,
                             std::span<const std::string> future_row_ids,
                             PredictivePolicy = {}) const;
  const char *method_id() const noexcept {
    return "proper-Gaussian-joint-predictive/retained-covariance/v1";
  }

private:
  void invalidate() noexcept;
  std::optional<Gaussian> predictive_;
  PredictiveMetadata metadata_;
  // Empty invalid owners must not allocate even metadata default labels.
  Metadata training_source_{"", {}, "", "", "", "", "", "", "", "", ""},
      future_noise_{"", {}, "", "", "", "", "", "", "", "", ""};
  std::string prior_identity_, training_design_identity_, parameter_measure_,
      prior_dependence_identity_;
  std::vector<std::string> shared_nuisance_ids_;
  std::vector<double> training_values_, response_, mean_, mean_error_,
      covariance_error_;
  PredictivePolicy policy_;
  long double inverse_norm_estimate_ = 0, covariance_perturbation_ = 0;
  double relative_covariance_error_ = 0;
  std::size_t preparation_work_units_ = 0;
  DensityStatus status_ = DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
};
} // namespace irred::statistics
