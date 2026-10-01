#pragma once
#include "irred/numerics.hpp"
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace irred::cosmology {
// Today's critical-density fractions. Radiation is massless; baryons are a
// subset of pressureless matter and photons a subset of radiation. Flatness
// fixes Lambda. No temperature, Neff, ionization or drag-epoch inference.
struct EarlyFlatModel {
  double h0_km_s_mpc, omega_m, omega_r, omega_b, omega_gamma;
};
struct SoundHorizonRequest {
  EarlyFlatModel model;
  double z_drag;
  std::string drag_origin;
};
struct SoundHorizonPolicy {
  double absolute_tolerance_mpc, relative_tolerance;
  std::size_t maximum_callbacks_per_point;
  unsigned maximum_depth;
  std::size_t maximum_points, maximum_total_callbacks, maximum_native_bytes;
};
struct SoundHorizonRow {
  SoundHorizonRequest source;
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> sound_horizon_mpc;
  // Scaled quadrature estimator plus integral-storage/final-cast diagnostics;
  // not a proven bound on quadrature or callback arithmetic.
  double error_estimate_mpc = 0;
  std::size_t callbacks = 0;
};
struct SoundHorizonBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::vector<SoundHorizonRow> rows;
  std::size_t callbacks = 0;
};
inline constexpr std::string_view sound_horizon_model_id =
    "flat-pressureless-matter-massless-radiation-lambda";
inline constexpr std::string_view sound_horizon_equation_id =
    "conditional-tight-coupling-sound-horizon-scale-factor";
// z_drag=0 is a mathematical control of this prescribed approximation, not a
// claim that actual late-time radiation and baryons remain tightly coupled.
// Conservative simultaneous owned payload for the supported STL, excluding
// borrowed requests, stack recursion, allocator bookkeeping and RSS. Empty
// batches have zero dynamic payload; the returned owner itself is stack storage.
// Evaluation requires FE_TONEAREST and long-double digits>=64/max_exponent>=16384.
std::optional<std::size_t> sound_horizon_payload_bound(
    std::size_t point_count, std::size_t origin_bytes_including_terminators) noexcept;
std::optional<std::size_t> sound_horizon_payload_bound(
    std::span<const SoundHorizonRequest>) noexcept;
SoundHorizonBatch evaluate_sound_horizon(
    std::span<const SoundHorizonRequest>, SoundHorizonPolicy);
} // namespace irred::cosmology
