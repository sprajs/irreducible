#pragma once
#include "irred/statistics.hpp"
namespace irred::statistics {
struct ParameterPrior {
  std::vector<std::string> ordered_parameter_ids, parameter_units,
      shared_nuisance_ids;
  std::vector<double> mean, covariance;
  std::string prior_identity, design_identity, residual_unit, parameter_measure,
      dependence_identity;
  bool noise_independence_declared = false;
};
struct PosteriorPolicy {
  std::size_t maximum_elements = 1000000;
  std::size_t maximum_payload_bytes = 256 * 1024 * 1024;
  // Conservative dimension-based preparation work admission, not elapsed time.
  std::size_t maximum_work_units = 100000000;
  double maximum_forward_sensitivity = 1e-10;
};
struct PosteriorMean {
  DensityStatus status = DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<double> value;
  std::vector<double> absolute_error_estimates;
  double backward_residual = 0, estimated_forward_sensitivity = 0;
};
// beta|r under a proper Gaussian prior, normalized in parameter_measure.
// Covariance/design/noise are fixed; this is distinct from estimator sampling.
class GaussianPosterior {
public:
  GaussianPosterior() = default;
  GaussianPosterior(const GaussianPosterior &) = delete;
  GaussianPosterior &operator=(const GaussianPosterior &) = delete;
  GaussianPosterior(GaussianPosterior &&) noexcept;
  GaussianPosterior &operator=(GaussianPosterior &&) noexcept;
  static GaussianPosterior prepare(Gaussian &&source,
                                   std::span<const double> design,
                                   std::span<const std::string> row_ids,
                                   ParameterPrior prior, PosteriorPolicy = {});
  static std::optional<std::size_t>
  preparation_payload_bound(const Gaussian &, const ParameterPrior &) noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  std::optional<std::size_t> evaluation_payload_bound() const noexcept;
  DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  const Gaussian &source() const noexcept { return source_; }
  const ParameterPrior &prior() const noexcept { return prior_; }
  std::span<const double> covariance() const noexcept {
    return posterior_.covariance();
  }
  std::span<const double> design() const noexcept { return design_; }
  double covariance_relative_error_estimate() const noexcept {
    return covariance_error_;
  }
  PosteriorMean condition(std::span<const double> residual,
                          std::span<const std::string> row_ids,
                          PosteriorPolicy = {}) const;
  GaussianResult log_density(std::span<const double> residual,
                             std::span<const std::string> row_ids,
                             std::span<const double> parameters,
                             std::span<const std::string> parameter_ids,
                             PosteriorPolicy = {}) const;
  const char *method_id() const noexcept {
    return "proper-Gaussian-parameter-posterior/whitened-precision/v1";
  }

private:
  Gaussian source_, posterior_;
  ParameterPrior prior_;
  std::vector<double> design_, precision_;
  std::vector<long double> whitened_design_;
  double covariance_error_ = 0, source_whitening_error_ = 0;
  long double source_inverse_norm_estimate_ = 0;
  DensityStatus status_ = DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
};
} // namespace irred::statistics
