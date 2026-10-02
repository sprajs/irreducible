#pragma once
#include "irred/detector.hpp"
namespace irred::detector {
enum class SelectionMeasure { joint_detection_record, selected_only };
struct Observation {
  double detection_threshold_adu = 0;
  bool detected = false;
  // With zero read noise use an integer electron count (discrete measure).
  // With positive read noise use measured_adu (continuous ADU measure).
  std::optional<std::uint32_t> electron_count;
  std::optional<double> measured_adu;
  SelectionMeasure measure = SelectionMeasure::joint_detection_record;
};
struct SelectionPolicy {
  std::size_t maximum_rows = 0, maximum_payload_bytes = 0;
  std::size_t maximum_poisson_terms = 256;
  double absolute_probability_allowance = 5e-13,
         relative_probability_allowance = 2e-11;
  // Continuous density is with respect to the caller's ADU coordinate.
  double absolute_density_allowance_per_adu = 5e-13,
         relative_density_allowance = 2e-11;
  // Dimensionless diagnostic for the returned logarithm.
  double absolute_log_allowance = 5e-13, scaled_log_allowance = 2e-11;
};
struct LikelihoodRow {
  Observation source;
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> log_value, detection_probability;
  double numerical_error_estimate = 0;
  // A structural zero has no finite log value, distinguished from failure.
  bool zero_probability = false;
};
struct LikelihoodBatch {
  numerics::Status status = numerics::Status::invalid_input;
  Input source;
  std::vector<LikelihoodRow> rows;
  std::size_t poisson_terms = 0;
  double omitted_poisson_tail_estimate = 0;
};
std::optional<std::size_t> likelihood_payload_bound(std::size_t) noexcept;
LikelihoodBatch likelihood(const Input &, std::span<const Observation>,
                           SelectionPolicy);
inline constexpr std::string_view selection_id =
    "DETECTOR/threshold-joint-or-selected-censoring/v1";
} // namespace irred::detector
