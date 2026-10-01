#pragma once
#include "irred/sound_horizon.hpp"
#include <array>
namespace irred::cosmology {
enum class EarlyLateOutput : unsigned {
  e,
  h_km_s_mpc,
  dh_mpc,
  dm_mpc,
  dl_mpc,
  dv_mpc,
  dm_over_rs,
  dh_over_rs,
  dv_over_rs,
  count
};
inline constexpr unsigned early_late_output_count =
    static_cast<unsigned>(EarlyLateOutput::count);
inline constexpr unsigned early_late_mask(EarlyLateOutput x) {
  const auto tag = static_cast<unsigned>(x);
  return tag < early_late_output_count ? 1u << tag : 0u;
}
struct EarlyLatePolicy {
  double absolute_tolerance_mpc, relative_tolerance;
  double absolute_tolerance_ratio, relative_tolerance_ratio;
  std::size_t maximum_callbacks_per_point, maximum_total_callbacks;
  unsigned maximum_depth;
  std::size_t maximum_points, maximum_native_bytes;
  SoundHorizonPolicy sound;
};
struct EarlyLateValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double error_estimate = 0;
};
struct EarlyLateRow {
  double redshift = 0;
  std::array<EarlyLateValue, early_late_output_count> outputs{};
  std::size_t callbacks = 0;
};
struct EarlyLateBatch {
  numerics::Status status = numerics::Status::invalid_input;
  SoundHorizonRequest source{};
  unsigned requested_outputs = 0;
  std::optional<SoundHorizonRow> ruler;
  std::vector<EarlyLateRow> rows;
  std::size_t callbacks = 0;
};
// One shared physical identity and supplied drag origin per batch. Unrequested
// outputs remain invalid_input and absent. Callback total includes the ruler.
std::optional<std::size_t>
early_late_payload_bound(std::size_t points, std::size_t origin_bytes) noexcept;
EarlyLateBatch evaluate_early_late(const SoundHorizonRequest &,
                                   std::span<const double> redshifts,
                                   unsigned requested_outputs, EarlyLatePolicy);
inline constexpr std::string_view early_late_equation_id =
    "flat-radiation-matter-lambda-log-redshift-distance";
} // namespace irred::cosmology
