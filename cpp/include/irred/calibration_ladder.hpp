#pragma once
#include "irred/gaussian_design.hpp"

namespace irred::calibration {
// Empirical supplied-shape model. Every input row is a synthetic control.
enum class RowKind {
  anchor_modulus,
  cepheid,
  calibrator_supernova,
  hubble_supernova,
  calibration_measurement
};
struct Row {
  RowKind kind = RowKind::anchor_modulus;
  std::string row_id, host_id, event_id;
  double calibration_response = 0;
  double log10_period_days = 0, metallicity_dex = 0;
  double reference_modulus_mag = 0;
};
struct Model {
  std::vector<std::string> ordered_host_ids;
  std::vector<Row> rows;
  double metallicity_reference_dex = 0;
  double h_reference_km_s_Mpc = 70;
  std::string distance_shape_identity, calibration_identity,
      dependence_identity;
  // Must describe a covariance conditional on delta, excluding a marginalized
  // contribution of the same shared calibration uncertainty.
  std::string conditional_covariance_identity;
};
struct Prediction {
  statistics::DensityStatus status = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<double> values_mag;
};
struct Recovery {
  statistics::DesignResult relative_fit;
  double h0_km_s_Mpc = 0; // usable only when relative_fit.status == finite
};
// Coefficients: ordered host mu, M_Cep, b, gamma, M_SN, eta, delta.
// IDs are generated from the model, not matched by position alone.
std::vector<std::string> parameter_ids(const Model &);
Prediction predict(const Model &, std::span<const double> coefficients,
                   std::span<const std::string> ordered_parameter_ids,
                   statistics::DesignPolicy = {});
class Ladder {
public:
  Ladder() = default;
  Ladder(const Ladder &) = delete;
  Ladder &operator=(const Ladder &) = delete;
  Ladder(Ladder &&) noexcept = default;
  Ladder &operator=(Ladder &&other) noexcept {
    if (this != &other) {
      model_ = std::move(other.model_);
      profile_ = std::move(other.profile_);
      preparation_status_ = other.preparation_status_;
      preparation_numerical_status_ = other.preparation_numerical_status_;
    }
    return *this;
  }
  // Successful preparation consumes covariance once. Failure preserves it.
  static Ladder prepare(statistics::Gaussian &&, Model,
                        statistics::DesignPolicy = {});
  statistics::DensityStatus status() const noexcept {
    return profile_.status() == statistics::DensityStatus::invalid_input
               ? preparation_status_
               : profile_.status();
  }
  numerics::Status numerical_status() const noexcept {
    return profile_.status() == statistics::DensityStatus::invalid_input
               ? preparation_numerical_status_
               : profile_.numerical_status();
  }
  statistics::DesignRank rank() const noexcept { return profile_.rank(); }
  const Model &model() const noexcept { return model_; }
  const statistics::DesignMetadata &design_metadata() const noexcept {
    return profile_.design_metadata();
  }
  const statistics::Metadata &observation_metadata() const noexcept {
    return profile_.metadata();
  }
  Recovery fit(std::span<const double> measured_mag,
               std::span<const std::string> ordered_row_ids,
               statistics::DesignPolicy = {}) const;
  Prediction predict(std::span<const double> coefficients,
                     std::span<const std::string> ordered_parameter_ids,
                     statistics::DesignPolicy = {}) const;

private:
  statistics::DensityStatus preparation_status_ =
      statistics::DensityStatus::invalid_input;
  numerics::Status preparation_numerical_status_ =
      numerics::Status::invalid_input;
  Model model_;
  statistics::DesignProfile profile_;
};
} // namespace irred::calibration
