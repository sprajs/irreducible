#pragma once
#include "irred/statistics.hpp"
namespace irred::statistics {
struct CalibrationPrior {
  std::vector<std::string> ordered_parameter_ids, parameter_units;
  std::vector<double> mean, covariance;
  std::string prior_identity, response_identity, residual_unit,
      calibration_identity, dependence_identity, measure_identity;
  bool noise_independence_declared = false;
};
struct CalibrationPolicy {
  std::size_t maximum_elements = 1000000;
  std::size_t maximum_payload_bytes = 256 * 1024 * 1024;
  double maximum_forward_sensitivity = 1e-10;
};
struct CalibrationResult : GaussianResult {
  double quadratic_rounding_estimate = 0, log_determinant_rounding_estimate = 0,
         log_density_rounding_estimate = 0;
};
// Normalized observed residual density; no posterior or model evidence.
class CorrelatedCalibration {
public:
  CorrelatedCalibration() = default;
  CorrelatedCalibration(const CorrelatedCalibration &) = delete;
  CorrelatedCalibration &operator=(const CorrelatedCalibration &) = delete;
  CorrelatedCalibration(CorrelatedCalibration &&) noexcept;
  CorrelatedCalibration &operator=(CorrelatedCalibration &&) noexcept;
  static CorrelatedCalibration prepare(Gaussian &&source,
                                       std::span<const double> response,
                                       std::span<const std::string> row_ids,
                                       CalibrationPrior prior,
                                       CalibrationPolicy policy = {});
  static std::optional<std::size_t>
  preparation_payload_bound(const Gaussian &,
                            const CalibrationPrior &) noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  std::optional<std::size_t> evaluation_payload_bound() const noexcept;
  DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  const Gaussian &source() const noexcept { return source_; }
  const CalibrationPrior &prior() const noexcept { return prior_; }
  std::span<const double> response() const noexcept { return response_; }
  std::span<const double> mean_shift() const noexcept { return shift_; }
  std::span<const double> effective_covariance() const noexcept {
    return effective_.covariance();
  }
  const char *method_id() const noexcept {
    return "correlated-proper-calibration/normalized-observed-density/v1";
  }
  numerics::Arithmetic arithmetic() const noexcept {
    return effective_.arithmetic();
  }
  CalibrationResult evaluate(std::span<const double> residual,
                             std::span<const std::string> row_ids,
                             CalibrationPolicy policy = {}) const;

private:
  double inverse_norm_estimate_ = 0, covariance_rounding_eta_ = 0;
  long double mean_rounding_inf_ = 0, mean_rounding_l1_ = 0;
  Gaussian source_, effective_;
  CalibrationPrior prior_;
  std::vector<double> response_, shift_;
  DensityStatus status_ = DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
};
} // namespace irred::statistics
