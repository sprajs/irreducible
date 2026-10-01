#pragma once
#include "irred/photometry.hpp"
namespace irred::photometry {
inline constexpr std::string_view sampled_model_id = "piecewise_linear_rest_luminosity_observed_optical_passband";
struct SampledSpectrum {
  std::span<const double> wavelength_metre, luminosity_watt_per_metre;
};
struct SampledPassband {
  std::span<const double> wavelength_metre, optical_transmission;
};
struct SampledInput {
  SampledSpectrum spectrum;
  SampledPassband passband;
  double luminosity_distance_metre = 0, redshift = 0;
  double collecting_area_square_metre = 0, observer_exposure_second = 0;
};
struct SampledPolicy {
  std::uint32_t requested_outputs = 0;
  // Total source + passband knots; conservative merged-segment count is knots-3.
  std::size_t maximum_samples = 0, maximum_segments = 0;
};
struct SampledResult {
  numerics::Status admission_status = numerics::Status::invalid_input;
  Outcome flux_watt_per_square_metre, energy_joule, expected_photons;
};
// Exact declared piecewise-linear models, zero outside finite supports.
// Incident flux integrates over passband SUPPORT without transmission weighting.
// Synchronous borrowed arrays, no input retention or allocation; owned scalar results.
// Hard total-knot cap 65536; policy/work-limit failures contain no output payload.
SampledResult evaluate_sampled(const SampledInput &, SampledPolicy) noexcept;
} // namespace irred::photometry
