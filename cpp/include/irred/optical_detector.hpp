#pragma once
#include "irred/detector_selection.hpp"
#include "irred/photometry_calibration.hpp"

namespace irred::photometry {
inline constexpr std::string_view optical_detector_model_id =
    "finite_joint_optical_state_conditional_independent_detector_all_band_selection_v1";
enum class ConditionalDetectorLaw { unspecified, independent_poisson_and_read };
// N and observer exposure come from the corresponding optical band/state.
struct DetectorBand {
  detector::PhotonLaw photon_law = detector::PhotonLaw::unspecified;
  double quantum_efficiency = 0, background_expected_electrons = 0;
  double dark_electrons_per_second = 0, read_noise_rms_electrons = 0;
  double gain_electrons_per_adu = 1, bias_adu = 0;
};
struct OpticalDetectorRecord {
  std::string_view id;
  std::span<const detector::Observation> bands;
  detector::SelectionMeasure measure = detector::SelectionMeasure::joint_detection_record;
};
struct OpticalDetectorInput {
  CalibrationInput optical;
  std::span<const DetectorBand> detector_bands; // exact optical band order
  std::span<const OpticalDetectorRecord> records;
  ConditionalDetectorLaw conditional_law = ConditionalDetectorLaw::unspecified;
  std::string_view detector_origin, observation_origin, conditional_law_origin;
};
struct OpticalDetectorPolicy {
  std::size_t maximum_states = 1024, maximum_bands = 64, maximum_records = 1024;
  std::size_t maximum_state_band_records = 65536;
  std::size_t maximum_total_optical_knots = 1048576;
  std::size_t maximum_optical_payload_bytes = 67108864;
  std::size_t maximum_payload_bytes = 134217728;
  detector::SelectionPolicy detector_policy{2048, 1048576};
  double absolute_log_allowance = 5e-12, scaled_log_allowance = 1e-8;
};
struct OpticalDetectorAttempt {
  std::string state_id, band_id;
  // Within each batch: original records, then threshold-only nondetection probes.
  detector::LikelihoodBatch nominal;
  std::optional<detector::LikelihoodBatch> lower_photons, upper_photons;
};
struct OpticalDetectorConditional {
  std::string state_id, record_id;
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> log_record, log_all_band_detection;
  bool zero_record = false, zero_all_band_detection = false;
  double record_log_error_estimate = 0, detection_log_error_estimate = 0;
};
struct OpticalDetectorMixture {
  std::string record_id;
  detector::SelectionMeasure measure = detector::SelectionMeasure::joint_detection_record;
  numerics::Status status = numerics::Status::invalid_input;
  // Joint mixture always retained, including when selected-only is requested.
  std::optional<double> log_joint_record, log_all_band_detection, log_value;
  bool zero_probability = false, zero_all_band_detection = false;
  double numerical_error_estimate = 0, joint_record_log_error_estimate = 0,
         detection_log_error_estimate = 0;
};
struct OpticalDetectorResult {
  numerics::Status status = numerics::Status::invalid_input;
  CalibrationResult optical;
  std::string detector_origin, observation_origin, conditional_law_origin;
  ConditionalDetectorLaw conditional_law = ConditionalDetectorLaw::unspecified;
  std::vector<OpticalDetectorAttempt> attempts; // state then band
  std::vector<OpticalDetectorConditional> conditional; // state then record
  std::vector<OpticalDetectorMixture> records; // original record order
  std::size_t poisson_terms = 0; // all nominal/endpoint attempts
};
// Synchronous borrowed inputs; result owns scalar sources, IDs and diagnostics.
// E=all bands detected. State masses and every required state refusal survive.
OpticalDetectorResult evaluate_optical_detector(const OpticalDetectorInput &,
                                                OpticalDetectorPolicy = {});
} // namespace irred::photometry
