#pragma once
#include "irred/statistics.hpp"
namespace irred::statistics {
// evaluate is an explicitly relative quadratic profile, never a density or
// evidence. GaussianBox separately owns a normalized finite-box completion.
enum class DesignRank {
  unassessed,
  deficient,
  unresolved,
  full_within_conditioning_contract
};
struct DesignMetadata {
  std::vector<std::string> ordered_parameter_ids, parameter_units,
      shared_nuisance_ids;
  std::string residual_unit, design_identity, dependence_identity;
};
struct DesignPolicy {
  std::size_t maximum_elements = 1000000;
  std::size_t maximum_payload_bytes = 256 * 1024 * 1024;
  double maximum_forward_sensitivity = 1e-10;
};
struct DesignResult {
  DensityStatus status = DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<double> coefficients, adjusted_residuals;
  double quadratic = 0, relative_log_score = 0;
  double normalized_normal_equation_residual = 0;
  double covariance_whitening_backward_residual = 0,
         covariance_whitening_rounding_estimate = 0;
};
// Declared weights map each named parameter unit into output_unit. These are
// caller provenance, not parsed conversions or a prior on fitted coefficients.
struct LinearFunctionalMetadata {
  std::vector<std::string> weight_units;
  std::string functional_identity, output_unit;
};
struct EstimatorVarianceResult {
  DensityStatus status = DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  double variance = 0, triangular_backward_residual = 0,
         estimated_forward_sensitivity = 0, output_rounding_error_relative = 0;
  LinearFunctionalMetadata metadata;
  const char *method_id = "retained-qr-linear-estimator-variance/v1";
  const char *sampling_law_id = "fixed-design-gaussian-observation-noise/v1";
};
class DesignProfile {
public:
  DesignProfile() = default;
  DesignProfile(const DesignProfile &) = delete;
  DesignProfile &operator=(const DesignProfile &) = delete;
  DesignProfile(DesignProfile &&) noexcept;
  DesignProfile &operator=(DesignProfile &&) noexcept;
  // X is row-major; p is the number of explicitly ordered parameters (>=2).
  // Failure leaves source usable; successful preparation consumes it.
  static DesignProfile prepare(Gaussian &&source, std::span<const double> x,
                               std::span<const std::string> ordered_row_ids,
                               DesignMetadata metadata,
                               DesignPolicy policy = {});
  static std::optional<std::size_t>
  preparation_payload_bound(const Gaussian &, std::size_t p,
                            const DesignMetadata &) noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  std::optional<std::size_t> evaluation_payload_bound() const noexcept;
  DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  DesignRank rank() const noexcept { return rank_; }
  double equilibrated_triangular_condition_inf() const noexcept {
    return triangular_condition_;
  }
  double equilibrated_transpose_triangular_condition_inf() const noexcept {
    return transpose_triangular_condition_;
  }
  const char *method_id() const noexcept {
    return "retained-whitened-pivoted-householder-qr/v1";
  }
  const char *qr_arithmetic_id() const noexcept { return "longdouble-cpu/v1"; }
  const Metadata &metadata() const noexcept { return gaussian_.metadata(); }
  const DesignMetadata &design_metadata() const noexcept {
    return design_metadata_;
  }
  const std::vector<SelectionRecord> &selection_history() const noexcept {
    return gaussian_.selection_history();
  }
  numerics::Arithmetic arithmetic() const noexcept {
    return gaussian_.arithmetic();
  }
  DesignResult evaluate(std::span<const double> residual,
                        std::span<const std::string> ordered_row_ids,
                        DesignPolicy policy = {}) const;
  // Conditional variance of the ideal linear estimator under supplied fixed
  // X,C and Gaussian observation noise. Not a posterior or rounding-noise law.
  EstimatorVarianceResult
  estimator_variance(std::span<const double> weights,
                     std::span<const std::string> ordered_parameter_ids,
                     LinearFunctionalMetadata, DesignPolicy = {}) const;
  std::optional<std::size_t> estimator_variance_payload_bound(
      const LinearFunctionalMetadata &) const noexcept;

private:
  friend class GaussianBox;
  Gaussian gaussian_;
  DesignMetadata design_metadata_;
  std::vector<double> x_;
  // A checks actual post-cast stationarity; QR contains R and reflector tails.
  std::vector<long double> whitened_design_, qr_, scales_, tau_;
  std::vector<std::size_t> pivot_;
  DensityStatus status_ = DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
  DesignRank rank_ = DesignRank::unassessed;
  double triangular_condition_ = 0, transpose_triangular_condition_ = 0,
         preparation_sensitivity_ = 0;
};
} // namespace irred::statistics
