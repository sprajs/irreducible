#pragma once
#include "irred/early_late.hpp"
namespace irred::lensing {
// Synthetic axial point source. Angles are radians, sigma_v is km/s, flux is
// an arbitrary fixed linear unit. For lambda!=1 sigma_v is the base SIS
// template parameter, not observed stellar kinematics; beta is the transformed
// source.
struct SISSource {
  cosmology::SoundHorizonRequest background;
  double lens_redshift, source_redshift, sigma_v_km_s, beta_radians;
  double source_flux, psf_sigma_radians, mass_sheet_lambda = 1;
};
struct SISPreparationPolicy {
  cosmology::EarlyLatePolicy geometry{
      1e-9,
      2e-11,
      1e-11,
      5e-11,
      1000000,
      4000000,
      40,
      2,
      16 * 1024 * 1024,
      {1e-9, 2e-11, 0, 40, 1, 0, 16 * 1024 * 1024}};
  std::size_t maximum_native_bytes = 16 * 1024 * 1024;
};
struct LensValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double absolute_error_estimate = 0;
};
inline constexpr unsigned lens_positions = 1, lens_magnifications = 2,
                          lens_fluxes = 4, lens_delays = 8;
struct SISImage {
  LensValue x_radians, signed_magnification, flux, relative_delay_seconds;
  int parity = 0;
};
struct SISPrediction {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested = 0;
  std::size_t image_count = 0;
  std::array<SISImage, 2> images;
};
struct PixelRectangle {
  double x_low, x_high, y_low, y_high;
};
struct SISPixelPolicy {
  double relative_tolerance = 1e-8;
  std::size_t maximum_pixels = 4096, maximum_native_bytes = 16 * 1024 * 1024,
              maximum_image_pixel_terms = 8192;
};
struct SISPixelRow {
  LensValue flux;
  std::size_t image_pixel_terms = 0, cdf_evaluations = 0;
};
struct SISPixelBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::vector<SISPixelRow> rows;
  std::size_t image_pixel_terms = 0, cdf_evaluations = 0;
};
class SISThinLens {
public:
  SISThinLens() = default;
  SISThinLens(const SISThinLens &) = default;
  SISThinLens &operator=(const SISThinLens &) = default;
  SISThinLens(SISThinLens &&) noexcept;
  SISThinLens &operator=(SISThinLens &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const SISSource *source() const noexcept {
    return source_ ? &*source_ : nullptr;
  }
  std::size_t preparation_callbacks() const noexcept { return callbacks_; }
  // Dd, Ds, Dds (physical Mpc), then base thetaE (radians).
  const std::array<LensValue, 4> &geometry() const noexcept {
    return geometry_;
  }
  SISPrediction predict(unsigned outputs,
                        double relative_tolerance = 1e-8) const;
  SISPixelBatch pixels(std::span<const PixelRectangle>,
                       SISPixelPolicy = {}) const;

private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<SISSource> source_;
  std::array<LensValue, 4> geometry_;
  std::size_t callbacks_ = 0, count_ = 0;
  std::array<long double, 2> x_{}, x_error_{}, mu_{}, mu_error_{}, flux_{},
      flux_error_{};
  long double delay_ = 0, delay_error_ = 0;
  friend SISThinLens prepare_sis_thin_lens(const SISSource &,
                                           SISPreparationPolicy);
};
SISThinLens prepare_sis_thin_lens(const SISSource &, SISPreparationPolicy = {});
std::optional<std::size_t> sis_payload_bound(std::size_t pixels,
                                             std::size_t origin_bytes) noexcept;
inline constexpr std::string_view sis_thin_lens_id =
    "synthetic-SIS-nonnegative-mass-sheet-point-source-Gaussian-rectangles/v1";
} // namespace irred::lensing
