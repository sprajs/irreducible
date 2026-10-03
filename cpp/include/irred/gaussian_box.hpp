#pragma once
#include "irred/gaussian_design.hpp"

namespace irred::statistics {
// Endpoints are exact supplied binary64 values in the original coordinate
// order. Fixed coordinates are point masses OUTSIDE this active measure.
struct BoxSupport {
  std::vector<std::string> ordered_parameter_ids;
  std::vector<double> lower, upper;
  std::string parameter_measure, prior_identity, fixed_coordinate_provenance;
};
struct BoxPolicy {
  DesignPolicy design;
  double maximum_log_probability_width = 1e-9;
  double maximum_quantile_width = 1e-8;
  // Standard-normal quadrature allocation only. Reported conditional-box CDF
  // diagnostics additionally include the excluded-mass bound.
  double cdf_absolute_radius = 2e-11;
  std::size_t maximum_cdf_nodes = 8192;
  std::size_t maximum_bisections = 96;
  // All CDF attempts, including unsuccessful refinement/inversion, count.
  std::size_t maximum_cdf_evaluations = 256;
};
struct BoxInterval { double lower = 0, upper = 0; };
struct BoxMarginalRequest {
  std::size_t active_parameter_index = 0;
  double cumulative_probability = .5;
};
enum class BoxStage {
  unassessed, gaussian_completion, endpoint_margins, rectangle_enclosure,
  normalization_enclosures, quantile_enclosure, complete
};
// Attempted completion gate, independent of the last available output stage.
// A refusal before completion still identifies the inherited numerical gate.
enum class BoxCompletionStep {
  unassessed, profile_evaluation, source_whitening, qr_completion,
  marginal_variances, covariance_determinant, complete
};
struct GaussianBoxResult {
  DensityStatus status = DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  BoxStage stage = BoxStage::unassessed;
  BoxCompletionStep completion_step = BoxCompletionStep::unassessed;
  // Present only while attempting a marginal variance in active source order.
  std::optional<std::size_t> completion_parameter_index;
  // The operation sets availability explicitly. Refusals retain earned
  // diagnostics, while zero-initialized absent fields are never predictions.
  bool gaussian_completion_available = false, endpoint_margins_available = false,
      rectangle_enclosure_available = false,
      normalization_enclosures_available = false,
      quantile_enclosure_available = false,
      endpoint_cdf_enclosures_available = false;
  // These completion inputs have empirical arithmetic diagnostics, not a
  // certified enclosure against the original supplied X,C,y.
  std::vector<double> unboxed_mean, unboxed_variance;
  double minimum_quadratic = 0, reported_profile_quadratic = 0;
  double log_design_precision_determinant = 0;
  double source_covariance_log_determinant = 0;
  double maximum_variance_sensitivity = 0, profile_stationarity = 0;
  // Conditional enclosures treat the reported Gaussian completion as exact.
  std::vector<BoxInterval> standardized_lower_margins,
      standardized_upper_margins;
  double excluded_mass_upper = 0;
  BoxInterval box_probability, log_box_probability, log_prior_volume;
  BoxInterval log_relative_box_integral, log_prior_normalized_relative_evidence,
      log_observation_normalized_evidence;
  BoxInterval requested_quantile, lower_endpoint_cdf, upper_endpoint_cdf;
  std::size_t cdf_node_evaluations = 0, cdf_evaluations = 0,
      bisections = 0;
  const char *method_id = "retained-qr-rational-tail-box-enclosure/v1";
  const char *enclosure_scope =
      "conditional-on-reported-gaussian-completion;not-original-input-certificate/v1";
};

// Minimal C++ consumer; no CLI/C ABI, Gaussian prior, constrained optimizer or
// chain sampler. Retains the original QR and covariance factor once.
class GaussianBox {
public:
  GaussianBox() = default;
  GaussianBox(const GaussianBox &) = delete;
  GaussianBox &operator=(const GaussianBox &) = delete;
  GaussianBox(GaussianBox &&) noexcept;
  GaussianBox &operator=(GaussianBox &&) noexcept;
  // Failure preserves the DesignProfile; success consumes it.
  static GaussianBox prepare(DesignProfile &&, BoxSupport, BoxPolicy = {});
  static std::optional<std::size_t> preparation_payload_bound(
      const DesignProfile &, const BoxSupport &) noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  std::optional<std::size_t> evaluation_payload_bound() const noexcept;
  DensityStatus status() const noexcept { return status_; }
  const BoxSupport &support() const noexcept { return support_; }
  GaussianBoxResult evaluate(std::span<const double> residual,
      std::span<const std::string> ordered_row_ids, BoxMarginalRequest,
      BoxPolicy = {}) const;
private:
  DesignProfile profile_;
  BoxSupport support_;
  DensityStatus status_ = DensityStatus::invalid_input;
};
} // namespace irred::statistics
