#pragma once
#include "irred/sampled_photometry.hpp"
#include <memory>
namespace irred::photometry {
inline constexpr std::string_view temporal_model_id =
    "finite_bilinear_rest_spectral_time_zero_outside_full_observer_exposure_"
    "mean";
inline constexpr std::string_view temporal_method_id =
    "linear_observer_time_bernstein_mean_then_shared_sampled_radiometry";
struct TemporalGrid {
  std::string_view id, source_origin;
  std::span<const double> rest_time_second, rest_wavelength_metre;
  // Rectangular time-major grid: [time_index * wavelength_count +
  // wavelength_index]. Exact declared linear functions in rest time and
  // wavelength, zero outside.
  std::span<const double> luminosity_watt_per_metre;
};
struct TemporalBand {
  std::string_view id, calibration_origin;
  SampledPassband passband;
};
struct TemporalPreparationPolicy {
  std::size_t maximum_grids = 4096, maximum_bands = 4096;
  std::size_t maximum_total_knots = 65536, maximum_grid_values = 4000000;
  std::size_t maximum_native_bytes = 1ull << 30;
};
struct TemporalExposure {
  std::size_t grid_index = 0, band_index = 0;
  double luminosity_distance_metre = 0, redshift = 0;
  double collecting_area_square_metre = 0;
  // Observer epoch of rest time zero, and absolute observer exposure interval.
  double source_epoch_observer_second = 0;
  double observer_lower_second = 0, observer_upper_second = 0;
};
struct TemporalPolicy {
  std::uint32_t requested_outputs = 7;
  std::size_t maximum_rows = 65536, maximum_native_bytes = 1ull << 30;
  // Cumulative conservative temporal cells + shared spectral segment admission.
  std::size_t maximum_segment_work = 4000000;
};
enum class TemporalCoverage { unassessed, no_overlap, partial, full };
struct TemporalRow {
  TemporalExposure source;
  numerics::Status admission_status = numerics::Status::invalid_input;
  TemporalCoverage coverage = TemporalCoverage::unassessed;
  // Wide diagnostics preserve positive coverage too small for binary64.
  long double observer_duration_second = 0, covered_observer_second = 0,
              covered_fraction = 0;
  std::size_t segment_work = 0;
  Outcome mean_flux_watt_per_square_metre, energy_joule, expected_photons;
};
struct TemporalBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::size_t segment_work = 0;
  std::vector<TemporalRow> rows;
};
struct TemporalStorage;
class PreparedTemporal {
public:
  PreparedTemporal() = default;
  PreparedTemporal(const PreparedTemporal &) = default;
  PreparedTemporal &operator=(const PreparedTemporal &) = default;
  PreparedTemporal(PreparedTemporal &&) noexcept;
  PreparedTemporal &operator=(PreparedTemporal &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  std::size_t grid_count() const noexcept;
  std::size_t band_count() const noexcept;
  std::optional<TemporalGrid> grid(std::size_t) const noexcept;
  std::optional<TemporalBand> band(std::size_t) const noexcept;
  // Immutable storage shared by copies; failed/moved owners have no payload.
  // Owner + native storage/requested vector payload; control-block bookkeeping,
  // borrowed caller spans, allocator overhead and RSS are excluded.
  std::optional<std::size_t> retained_payload_bytes() const noexcept;

private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::shared_ptr<const TemporalStorage> storage_;
  friend PreparedTemporal prepare_temporal(std::span<const TemporalGrid>,
                                           std::span<const TemporalBand>,
                                           TemporalPreparationPolicy);
  friend TemporalBatch evaluate_temporal(const PreparedTemporal &,
                                         std::span<const TemporalExposure>,
                                         TemporalPolicy);
};
// Bounds before scans/copying. Admitted axes have >=2 finite strictly
// increasing samples; wavelength>0, grid L>=0, observed optical T in [0,1].
// Calibration is supplied/fixed. Invalid preparation returns no source payload.
// Allocation failures are ordinary C++ exceptions; there is no C ABI or
// exception crossing.
PreparedTemporal prepare_temporal(std::span<const TemporalGrid>,
                                  std::span<const TemporalBand>,
                                  TemporalPreparationPolicy = {});
// Full observer-interval mean, not a mean conditional on covered time. Observer
// duration must be positive and representable with admitted relative rounding.
// D_L>0,z>=0,A>=0, all scalars finite. Partial/no time support is explicit.
// Row-level failures preserve other rows and independent requested groups.
// Coarse batch, one reused bounded scratch, no per-segment sampled calls.
TemporalBatch evaluate_temporal(const PreparedTemporal &,
                                std::span<const TemporalExposure>,
                                TemporalPolicy = {});
} // namespace irred::photometry
