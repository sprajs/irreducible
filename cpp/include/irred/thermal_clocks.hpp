#pragma once
#include "irred/thermal_neutrino.hpp"
#include <array>
namespace irred::cosmology {
enum class ThermalClockOutput : unsigned {
  age_gyr, lookback_gyr, comoving_particle_horizon_mpc,
  proper_particle_horizon_mpc, count
};
inline constexpr unsigned thermal_clock_output_count=
    static_cast<unsigned>(ThermalClockOutput::count);
constexpr unsigned thermal_clock_mask(ThermalClockOutput output) {
  const auto tag=static_cast<unsigned>(output);
  return tag<thermal_clock_output_count ? 1u<<tag : 0;
}
struct ThermalClockPolicy {
  double absolute_tolerance_gyr=1e-8, absolute_tolerance_mpc=1e-8;
  double relative_tolerance=2e-10;
  std::size_t maximum_callbacks_per_point=20000000;
  std::size_t maximum_total_callbacks=100000000;
  std::size_t maximum_points=4096, maximum_native_bytes=16*1024*1024;
  unsigned maximum_depth=30, maximum_tail_refinements=80;
  ThermalPolicy thermal=[] { auto p=ThermalPolicy{};
    p.maximum_total_callbacks=100000000; return p; }();
};
struct ThermalClockTail {
  // Dimensionless H0*t or H0*eta contribution from the retained full [0,ae]
  // support. These empirical endpoint brackets are not certified bounds.
  long double end_scale_factor=0, midpoint=0, half_width=0;
  unsigned refinements=0;
};
struct ThermalClockRow {
  double scale_factor=0;
  std::array<ThermalBackgroundValue,thermal_clock_output_count> outputs{};
  std::optional<ThermalClockTail> age_tail, conformal_tail;
  std::size_t callbacks=0, outer_callbacks=0, momentum_callbacks=0;
};
struct ThermalClockBatch {
  numerics::Status status=numerics::Status::invalid_input;
  unsigned requested_outputs=0;
  std::vector<ThermalClockRow> rows;
  std::size_t callbacks=0, outer_callbacks=0, momentum_callbacks=0;
};
class ThermalClocks {
public:
  ThermalClocks()=default;
  ThermalClocks(const ThermalClocks &)=default;
  ThermalClocks &operator=(const ThermalClocks &)=default;
  ThermalClocks(ThermalClocks &&) noexcept;
  ThermalClocks &operator=(ThermalClocks &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const ThermalBackground &background() const noexcept { return background_; }
  ThermalClockBatch evaluate(std::span<const double> scale_factors,
                            unsigned requested_outputs,
                            ThermalClockPolicy={}) const;
private:
  numerics::Status status_=numerics::Status::invalid_input;
  ThermalBackground background_;
  ThermalStressSourceEstimate source_estimate_;
  long double radiation_=0, radiation_error_=0;
  long double lambda_=0, omega_species_today_=0;
  friend ThermalClocks prepare_thermal_clocks(const ThermalBackground &,
      ThermalClockPolicy,ThermalStressSourceEstimate);
};
ThermalClocks prepare_thermal_clocks(const ThermalBackground &,
    ThermalClockPolicy={},ThermalStressSourceEstimate={});
std::optional<std::size_t> thermal_clocks_payload_bound(
    std::size_t points,std::size_t species) noexcept;
inline constexpr std::string_view thermal_clocks_model_id=
    "flat-thermal-FD-effective-extrapolation-proper-clock-particle-horizon/v1";
inline constexpr std::string_view thermal_clocks_method_id=
    "monotone-retained-early-tail-direct-lookback-scaled-adaptive-clock/v1";
inline constexpr std::string_view thermal_clocks_units_id=
    "Gyr-exact-Julian-year-SI-IAU-Mpc-horizon/v1";
} // namespace irred::cosmology
