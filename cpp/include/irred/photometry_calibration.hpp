#pragma once
#include "irred/sampled_photometry.hpp"
#include <string>
namespace irred::photometry {
inline constexpr std::string_view calibration_model_id =
    "finite_joint_optical_passband_calibration_v1";
struct CalibrationBand {
  std::string_view id;
  std::span<const double> observed_wavelength_metre;
  double collecting_area_square_metre = 0, observer_exposure_second = 0;
};
struct CalibrationState {
  std::string_view id;
  double relative_mass = 0;
  // One transmission array per band, in exactly the declared band order.
  std::span<const std::span<const double>> optical_transmission;
};
struct CalibrationInput {
  SampledSpectrum spectrum;
  double luminosity_distance_metre = 0, redshift = 0;
  std::span<const CalibrationBand> bands;
  std::span<const CalibrationState> states;
  std::string_view source_origin, calibration_origin, distribution_origin,
      dependence_origin;
};
struct CalibrationPolicy {
  std::uint32_t requested_outputs = collected_energy | transmitted_photons;
  std::size_t maximum_states = 1024, maximum_bands = 64;
  std::size_t maximum_total_knots = 1048576, maximum_output_bytes = 67108864;
  double mean_relative_sensitivity = 1e-10,
         covariance_relative_sensitivity = 1e-8;
};
struct CalibrationAxis {
  std::string band_id;
  std::uint32_t output;
}; // energy J; photons count
struct CalibrationAttempt {
  std::string state_id, band_id;
  SampledResult result;
};
struct CalibrationMoments {
  std::vector<double> mean,
      covariance; // axis order; covariance row major, product units
  std::vector<double> mean_empirical_sensitivity,
      covariance_empirical_sensitivity;
};
struct CalibrationResult {
  numerics::Status status = numerics::Status::invalid_input;
  std::string source_origin, calibration_origin, distribution_origin,
      dependence_origin;
  std::vector<CalibrationAxis> axes; // band order, then energy before photons
  std::vector<CalibrationAttempt>
      attempts; // state order, then band order; failures preserved
  std::optional<CalibrationMoments>
      moments; // absent if any required state/output fails
};
// Synchronous borrowed input. All masses are positive finite normal binary64;
// normalize once in wide arithmetic. Full joint law, never shot/detection
// noise. Explicit input equality or zero area/exposure/source witnesses can
// establish constant axes; equal rounded outputs alone cannot. No jitter or
// dropped states.
CalibrationResult evaluate_calibration(const CalibrationInput &,
                                       CalibrationPolicy = {});
} // namespace irred::photometry
