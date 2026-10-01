#pragma once
#include "irred/numerics.hpp"
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>
namespace irred::photometry {
inline constexpr std::string_view model_id = "constant_rest_luminosity_rectangular_band";
inline constexpr std::string_view constants_id = "si_2019_radiometric_definitions";
inline constexpr std::uint32_t incident_flux = 1, collected_energy = 2, transmitted_photons = 4;
enum class Availability : std::uint32_t { omitted, available, failed };
struct Input {
  double luminosity_watt_per_metre = 0;
  double rest_lower_metre = 0, rest_upper_metre = 0;
  double observed_lower_metre = 0, observed_upper_metre = 0;
  double luminosity_distance_metre = 0, redshift = 0;
  double collecting_area_square_metre = 0, optical_transmission = 0;
  double observer_exposure_second = 0;
};
struct Outcome {
  Availability availability = Availability::omitted;
  numerics::Status numerical_status = numerics::Status::ok;
  std::optional<double> value;
};
struct Row {
  Input source;
  numerics::Status admission_status = numerics::Status::invalid_input;
  Outcome flux_watt_per_square_metre, energy_joule, expected_photons;
};
struct Policy {
  std::uint32_t requested_outputs = 0;
  std::size_t maximum_rows = 0, maximum_output_bytes = 0;
};
struct Batch {
  numerics::Status status = numerics::Status::invalid_input;
  std::vector<Row> rows;
};
// Owns attempts/results; input buffers are borrowed only during the call.
// Bound covers successful Batch/Row allocation payload; no numerical scratch allocation.
// The ABI may allocate a fixed bounded empty work-limit diagnostic owner even at byte quota zero.
// Caller storage, allocator overhead and RSS are excluded.
std::optional<std::size_t> output_payload_bound(std::size_t rows) noexcept;
Batch evaluate(std::span<const Input>, Policy);
} // namespace irred::photometry
