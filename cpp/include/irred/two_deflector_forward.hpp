#pragma once
#include "irred/numerics.hpp"
#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>

namespace irred::lensing {
struct ForwardPoint2 { double x = 0, y = 0; };
// Dimensionless u=theta/theta_scale, one lens plane. This is an elliptical
// POTENTIAL, not an elliptical-density SPLE or observed stellar dispersion.
struct SoftenedPotentialComponent {
  ForwardPoint2 center;
  double strength = 0, core = 0, axis_ratio = 1, angle_radians = 0;
};
struct ExtendedGaussianBrightness {
  ForwardPoint2 center;
  double major_width = 0, minor_width = 0, angle_radians = 0;
  double peak_electrons_per_second_per_radian_squared = 0;
};
struct TwoDeflectorScene {
  std::array<SoftenedPotentialComponent, 2> deflectors;
  ExtendedGaussianBrightness source;
  double shear_1 = 0, shear_2 = 0, mass_scale_lambda = 1;
  double theta_scale_radians = 0, exposure_seconds = 0;
  double uniform_sky_electrons_per_second_per_radian_squared = 0;
  // Supplied steady observer-band electron-brightness/response provenance.
  // This does not qualify an ACS exposure, arrival law or measured noise.
  std::string origin;
};
struct ShiftPSFComponent { ForwardPoint2 shift; double probability = 0; };
struct AffineDetectorCutout {
  ForwardPoint2 origin;
  // Row-major detector coordinates -> u. No WCS interpretation.
  std::array<double, 4> matrix{1, 0, 0, 1};
};
struct ForwardPixelRectangle {
  double x_low = 0, x_high = 0, y_low = 0, y_high = 0;
  std::uint64_t original_index = 0;
};
struct ForwardWork {
  std::size_t preparation_started = 0, pixels_started = 0;
  std::size_t index_comparisons_started = 0, psf_preflight_started = 0;
  // Includes preflight maps and every attempted PSF field term. Refused
  // next terms do not count as begun; failed_starts records that refusal.
  std::size_t field_samples_started = 0, deflector_evaluations_started = 0;
  std::size_t outer_integrations_started = 0, inner_integrations_started = 0;
  std::size_t outer_callbacks_started = 0, inner_callbacks_started = 0;
  std::size_t failed_starts = 0;
};
// All contributors are empirical diagnostics in electrons. In particular,
// the maximum sampled inner error is not a uniform unsampled error bound.
struct ForwardErrorAllocation {
  double inner_quadrature_electrons = 0, outer_quadrature_electrons = 0;
  double field_psf_arithmetic_electrons = 0, projection_area_electrons = 0;
  double total_electrons = 0;
};
struct TwoDeflectorPreparationPolicy {
  std::size_t maximum_origin_bytes = 65536;
  std::size_t maximum_native_bytes = 4 * 1024 * 1024;
};
struct ForwardPixelPolicy {
  // Positive absolute floor is required; relative tolerance is nonnegative.
  double absolute_tolerance_electrons = 1e-10, relative_tolerance = 1e-10;
  std::size_t maximum_pixels = 4096;
  std::size_t maximum_field_samples = 8000000;
  std::size_t maximum_field_samples_per_pixel = 500000;
  std::size_t maximum_native_bytes = 4 * 1024 * 1024;
  unsigned maximum_depth = 20;
};
struct ForwardPixelMean {
  std::uint64_t original_index = 0;
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> electrons;
  // Absent until all contributors are computed; a budget refusal retains the
  // computed record without exposing a mean.
  std::optional<ForwardErrorAllocation> errors;
  ForwardWork work;
};
struct ForwardPixelMeans {
  ForwardPixelMeans() = default;
  ForwardPixelMeans(const ForwardPixelMeans &) = delete;
  ForwardPixelMeans &operator=(const ForwardPixelMeans &) = delete;
  ForwardPixelMeans(ForwardPixelMeans &&) noexcept;
  ForwardPixelMeans &operator=(ForwardPixelMeans &&) noexcept;
  numerics::Status status = numerics::Status::invalid_input;
  std::span<const ForwardPixelMean> rows() const noexcept {
    return {rows_.get(), row_count_};
  }
  ForwardWork work;
  // Requested simultaneous owned payload, excluding borrowed input,
  // recursion stack, allocator bookkeeping and process RSS.
  std::optional<std::size_t> admitted_payload_bytes;
private:
  std::unique_ptr<ForwardPixelMean[]> rows_;
  std::size_t row_count_ = 0;
  friend class PreparedTwoDeflectorForward;
};
class PreparedTwoDeflectorForward {
public:
  PreparedTwoDeflectorForward() = default;
  PreparedTwoDeflectorForward(const PreparedTwoDeflectorForward &) = delete;
  PreparedTwoDeflectorForward &operator=(const PreparedTwoDeflectorForward &) = delete;
  PreparedTwoDeflectorForward(PreparedTwoDeflectorForward &&) noexcept;
  PreparedTwoDeflectorForward &operator=(PreparedTwoDeflectorForward &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const ForwardWork &preparation_work() const noexcept { return work_; }
  const TwoDeflectorScene *scene() const noexcept {
    return scene_ ? &*scene_ : nullptr;
  }
  ForwardPixelMeans means(std::span<const ForwardPixelRectangle>,
                          ForwardPixelPolicy = {}) const;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
private:
  void swap(PreparedTwoDeflectorForward &) noexcept;
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<TwoDeflectorScene> scene_;
  AffineDetectorCutout cutout_;
  std::array<ShiftPSFComponent, 9> psf_{};
  std::size_t psf_count_ = 0;
  // Symmetric row-major Q0,Q1,P prepared once.
  std::array<std::array<long double, 4>, 3> derived_{};
  long double determinant_ = 0;
  long double determinant_error_ = 0;
  ForwardWork work_;
  friend PreparedTwoDeflectorForward prepare_two_deflector_forward(
      TwoDeflectorScene, const AffineDetectorCutout &,
      std::span<const ShiftPSFComponent>, TwoDeflectorPreparationPolicy);
};
// Consume an acquired scene. Passing an lvalue copies before native admission;
// that caller-owned copy is not part of the native requested-payload contract.
PreparedTwoDeflectorForward prepare_two_deflector_forward(
    TwoDeflectorScene, const AffineDetectorCutout &,
    std::span<const ShiftPSFComponent>, TwoDeflectorPreparationPolicy = {});
// Source-derived conservative request charge, using actual origin CAPACITY.
std::optional<std::size_t> two_deflector_payload_bound(
    std::size_t pixels, std::size_t origin_capacity) noexcept;
inline constexpr std::string_view two_deflector_forward_id =
    "synthetic-single-plane-two-softened-elliptical-potentials-extended-Gaussian-finite-shift-PSF/v1";
inline constexpr std::string_view two_deflector_forward_method_id =
    "resolved-pixel-nested-adaptive-Simpson/empirical-diagnostics/v1";
inline constexpr std::string_view two_deflector_forward_arithmetic_id =
    "binary64-callback/longdouble-scene/FE_TONEAREST/v1";
} // namespace irred::lensing
