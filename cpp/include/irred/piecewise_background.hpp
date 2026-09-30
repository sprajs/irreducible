#pragma once
#include "irred/background.hpp"
#include <array>
#include <optional>
namespace irred::cosmology {
namespace detail {
struct PiecewiseRadialAccess;
}
inline constexpr std::array<double, 6> piecewise_q_edges{0, .1, .3, .6, 1, 2.5};
struct PiecewiseQParameters {
  double h0_km_s_mpc;
  std::array<double, 5> q;
  PiecewiseQParameters(double h0, std::array<double, 5> values)
      : h0_km_s_mpc(h0), q(values) {}
};
// These derivative statements apply to the compiled piecewise model only.
enum class QConvention : std::uint32_t {
  not_assessed,
  interior_constant_bin,
  right_limit_at_internal_jump,
  right_limit_at_zero,
  left_limit_at_final_endpoint
};
enum class JerkAvailability : std::uint32_t {
  not_assessed,
  ordinary_within_bin,
  one_sided_endpoint,
  unavailable_at_jump
};
struct PiecewiseGeometry {
  double expansion_E = 0, h_km_s_mpc = 0, radial_integral = 0;
  double radial_mpc = 0, transverse_mpc = 0, angular_diameter_mpc = 0,
         luminosity_mpc = 0;
  double dimensionless_luminosity_shape = 0, lookback_seconds = 0;
  double volume_mpc3_per_sr_per_redshift = 0;
};
struct PiecewiseSlot {
  Query source{};
  Status status = Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  PiecewiseGeometry geometry;
  std::optional<double> assigned_q, jerk;
  QConvention q_convention = QConvention::not_assessed;
  JerkAvailability jerk_availability = JerkAvailability::not_assessed;
  std::size_t bin = 0, segments_processed = 0;
  // Only at z0: same retained first bin, imposed over interval[0,.1).
  std::optional<double> q0_within_piecewise_model;
};
struct PiecewisePolicy {
  std::size_t maximum_queries = 4096;
  // Analytic admitted segment visits, not quadrature callback counts.
  std::size_t maximum_segment_visits = 20480;
};
struct PiecewiseBatch {
  Status status = Status::invalid_input;
  std::vector<PiecewiseSlot> slots;
  std::size_t segments_processed = 0;
};
class PiecewiseBackground {
public:
  Status status() const noexcept { return status_; }
  const PiecewiseQParameters &parameters() const noexcept {
    return parameters_;
  }
  static constexpr std::string_view model_id =
      "P01/fixed-five-bin-q-flat-kinematic/v1";
  static constexpr std::string_view constants_id = irred::constant_set_id;
  PiecewiseBatch evaluate_batch(std::span<const Query>, PiecewisePolicy) const;

private:
  PiecewiseQParameters parameters_{0, {}};
  Status status_ = Status::invalid_input;
  double hubble_distance_mpc_ = 0, hubble_time_seconds_ = 0;
  std::array<long double, 5> start_E_{};
  friend struct detail::PiecewiseRadialAccess;
  friend PiecewiseBackground prepare_piecewise_q(PiecewiseQParameters);
};
// Separate analytic provider; no legacy v1/v2 model-tag admission.
// No extrapolation, smoothing prior, age/rd or sampler.
// Finite qi[-3,2], expansion z[0,2.5], finite positive representable
// H0. Geometry's nonzero outputs must be normal binary64; exact z0 zeros
// allowed. Retained qi may be subnormal or zero: they are explicit
// dimensionless inputs, not physical-length outputs. Undefined jerk is absent,
// never NaN or finite0.
PiecewiseBackground prepare_piecewise_q(PiecewiseQParameters);
} // namespace irred::cosmology
