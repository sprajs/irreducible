#pragma once
#include "irred/gr_growth.hpp"
#include <string>
namespace irred::cosmology {
enum class Sigma8Convention : std::uint32_t {
  linear_pressureless_total_matter_top_hat_8_over_h_mpc,
  unknown = UINT32_MAX
};
enum class AmplitudeTreatment : std::uint32_t {
  fixed_supplied,
  unknown = UINT32_MAX
};
struct FixedSigma8Source {
  double sigma8_ref = 0, reference_scale_factor = 0;
  Sigma8Convention convention = Sigma8Convention::unknown;
  AmplitudeTreatment treatment = AmplitudeTreatment::unknown;
  std::string amplitude_identity, amplitude_provenance;
};
inline constexpr unsigned amplitude_sigma8 = 1, amplitude_f_sigma8 = 2;
struct GrowthAmplitudePolicy {
  GrowthPolicy growth;
  double absolute_tolerance = 1e-12, relative_tolerance = 5e-10;
  size_t maximum_string_bytes = 1024 * 1024,
         maximum_native_bytes = 16 * 1024 * 1024;
};
struct GrowthAmplitudeValue {
  Availability availability = Availability::not_requested;
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double absolute_error_estimate = 0;
};
struct GrowthAmplitudeRow {
  double scale_factor = 0;
  GrowthAmplitudeValue sigma8, f_sigma8;
  size_t callbacks = 0;
};
struct GrowthAmplitudeBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested = 0;
  std::vector<GrowthAmplitudeRow> rows;
  std::optional<GrowthRow> reference;
  size_t callbacks = 0;
};
numerics::Status fixed_sigma8_status(const FixedSigma8Source &) noexcept;
// Returned rows plus simultaneous serial reference/query provider scratch.
// Excludes borrowed source/queries, retained growth, allocator overhead and
// RSS.
std::optional<size_t> growth_amplitude_payload_bound(size_t points) noexcept;
GrowthAmplitudeBatch evaluate_growth_amplitude(const GRGrowth &,
                                               const FixedSigma8Source &,
                                               std::span<const double>,
                                               unsigned requested,
                                               GrowthAmplitudePolicy = {});
inline constexpr std::string_view growth_amplitude_id =
    "GR/flat-pressureless-supplied-linear-sigma8-reference/v1";
} // namespace irred::cosmology
