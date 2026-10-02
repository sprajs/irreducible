#pragma once
#include "irred/calibration_ladder.hpp"
#include "irred/gaussian_predictive.hpp"
namespace irred::calibration {
struct H0PosteriorProjection {
  statistics::DensityStatus status = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  double eta_mean = 0, eta_variance = 0, median_h0_km_s_Mpc = 0,
         expectation_h0_km_s_Mpc = 0, log_standard_deviation = 0;
  double eta_mean_absolute_error_estimate = 0,
         eta_variance_absolute_error_estimate = 0,
         median_relative_error_estimate = 0,
         expectation_relative_error_estimate = 0,
         log_standard_deviation_relative_error_estimate = 0;
};
// Proper-prior synthetic ladder. Observation values retain their offsets;
// offset subtraction must be exactly representable, otherwise evaluation fails.
class LadderPosterior {
public:
  LadderPosterior() = default;
  LadderPosterior(const LadderPosterior &) = delete;
  LadderPosterior &operator=(const LadderPosterior &) = delete;
  LadderPosterior(LadderPosterior &&) noexcept;
  LadderPosterior &operator=(LadderPosterior &&) noexcept;
  static LadderPosterior prepare(statistics::Gaussian &&, Model,
                                 statistics::ParameterPrior,
                                 statistics::PosteriorPolicy = {});
  static std::optional<std::size_t>
  preparation_payload_bound(const statistics::Gaussian &, const Model &,
                            const statistics::ParameterPrior &) noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  std::optional<std::size_t> evaluation_payload_bound() const noexcept;
  statistics::DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  const Model &model() const noexcept { return model_; }
  const Linearization &linearization() const noexcept { return linear_; }
  // Requires status()==finite; empty owners allocate no fallback parent.
  const statistics::GaussianPosterior &posterior() const noexcept {
    return *posterior_;
  }
  statistics::PosteriorMean condition(std::span<const double>,
                                      std::span<const std::string>,
                                      statistics::PosteriorPolicy = {}) const;
  statistics::GaussianResult
      log_density(std::span<const double>, std::span<const std::string>,
                  std::span<const double>, std::span<const std::string>,
                  statistics::PosteriorPolicy = {}) const;
  H0PosteriorProjection h0_projection(std::span<const double>,
                                      std::span<const std::string>,
                                      statistics::PosteriorPolicy = {}) const;

private:
  Model model_;
  Linearization linear_;
  std::optional<statistics::GaussianPosterior> posterior_;
  statistics::PosteriorPolicy policy_;
  statistics::DensityStatus status_ = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
};
class LadderPredictive {
public:
  LadderPredictive() = default;
  LadderPredictive(const LadderPredictive &) = delete;
  LadderPredictive &operator=(const LadderPredictive &) = delete;
  LadderPredictive(LadderPredictive &&) noexcept;
  LadderPredictive &operator=(LadderPredictive &&) noexcept;
  static LadderPredictive
  prepare(const LadderPosterior &, std::span<const double> training_mag,
          std::span<const std::string>,
          const statistics::Gaussian &future_noise, Model future_model,
          statistics::PredictiveMetadata, statistics::PredictivePolicy = {});
  static std::optional<std::size_t>
  preparation_payload_bound(const LadderPosterior &,
                            const statistics::Gaussian &, const Model &,
                            const statistics::PredictiveMetadata &) noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  std::optional<std::size_t> evaluation_payload_bound() const noexcept;
  statistics::DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  const Model &model() const noexcept { return model_; }
  const Linearization &linearization() const noexcept { return linear_; }
  const statistics::GaussianPredictive &predictive() const noexcept {
    return predictive_;
  }
  std::span<const double> mean() const noexcept { return mean_; }
  std::span<const double> mean_absolute_error_estimates() const noexcept {
    return mean_errors_;
  }
  std::span<const double> covariance() const noexcept {
    return predictive_.covariance();
  }
  statistics::GaussianResult
      log_density(std::span<const double>, std::span<const std::string>,
                  statistics::PredictivePolicy = {}) const;

private:
  Model model_;
  Linearization linear_;
  statistics::GaussianPredictive predictive_;
  std::vector<double> mean_, mean_errors_;
  statistics::PredictivePolicy policy_;
  statistics::DensityStatus status_ = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
};
} // namespace irred::calibration
