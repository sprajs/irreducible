#include "irred/sis_thin_lens.hpp"
#include "irred/quantities.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <numbers>
namespace irred::lensing {
namespace {
using W = long double;
using S = numerics::Status;
constexpr W eps = std::numeric_limits<double>::epsilon();
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
}
bool tolerance(double r) { return std::isfinite(r) && r > 0; }
void project(LensValue &out, W value, W error, double relative) {
  double v = (double)value;
  W e = error + std::abs(value - v) + 64 * eps * std::abs(value);
  if (!std::isfinite(value) || !std::isfinite(error) || error < 0 ||
      (value != 0 && !std::isnormal(v)) || !std::isfinite(v)) {
    out.status = S::outside_domain;
    return;
  }
  if (e > relative * std::abs(value) || !std::isfinite((double)e)) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out = {S::ok, v, (double)e};
}
// The tail branch avoids subtracting two CDF values near one. Error is a
// declared arithmetic diagnostic of platform erfc/erf, not a certified bound.
W normal_interval(W lo, W hi, W &error) {
  W a = lo / std::sqrt(2.L), b = hi / std::sqrt(2.L), p = 0,
    arithmetic_scale = 0;
  if (a >= 0) {
    W first = std::erfc(a), second = std::erfc(b);
    p = (first - second) / 2;
    arithmetic_scale = (first + second) / 2;
  } else if (b <= 0) {
    W first = std::erfc(-b), second = std::erfc(-a);
    p = (first - second) / 2;
    arithmetic_scale = (first + second) / 2;
  } else {
    W first = std::erf(b), second = std::erf(a);
    p = (first - second) / 2;
    arithmetic_scale = (std::abs(first) + std::abs(second)) / 2;
  }
  error =
      128 * std::numeric_limits<W>::epsilon() * arithmetic_scale + 64 * eps * p;
  return p;
}
} // namespace
std::optional<size_t> sis_payload_bound(size_t pixels,
                                        size_t origins) noexcept {
  irred::detail::PayloadAccounting b(
      sizeof(SISThinLens) + sizeof(SISPixelBatch) +
      4 * sizeof(cosmology::EarlyLateRow) +
      4 * sizeof(cosmology::SoundHorizonRequest) + 2048);
  b.add(pixels, sizeof(SISPixelRow));
  b.add(origins, 4);
  return b.result();
}
SISThinLens::SISThinLens(SISThinLens &&o) noexcept { *this = std::move(o); }
SISThinLens &SISThinLens::operator=(SISThinLens &&o) noexcept {
  if (this != &o) {
    status_ = o.status_;
    source_ = std::move(o.source_);
    geometry_ = o.geometry_;
    callbacks_ = o.callbacks_;
    count_ = o.count_;
    x_ = o.x_;
    x_error_ = o.x_error_;
    mu_ = o.mu_;
    mu_error_ = o.mu_error_;
    flux_ = o.flux_;
    flux_error_ = o.flux_error_;
    delay_ = o.delay_;
    delay_error_ = o.delay_error_;
    o.source_.reset();
    o.status_ = S::invalid_input;
    o.callbacks_ = o.count_ = 0;
    o.geometry_ = {};
  }
  return *this;
}
SISThinLens prepare_sis_thin_lens(const SISSource &source,
                                  SISPreparationPolicy p) {
  SISThinLens out;
  if (!arithmetic())
    return out;
  const auto origin = source.background.drag_origin.size();
  if (origin > 65536) {
    out.status_ = S::work_limit;
    return out;
  }
  auto bytes = sis_payload_bound(0, origin + 1);
  if (!bytes || *bytes > p.maximum_native_bytes ||
      *bytes > size_t(1024) * 1024 * 1024) {
    out.status_ = S::work_limit;
    return out;
  }
  for (double v : {source.lens_redshift, source.source_redshift,
                   source.sigma_v_km_s, source.beta_radians, source.source_flux,
                   source.psf_sigma_radians, source.mass_sheet_lambda})
    if (!std::isfinite(v)) {
      out.status_ = S::nonfinite_input;
      return out;
    }
  W zl = source.lens_redshift, zs = source.source_redshift,
    sigma = source.sigma_v_km_s, beta = source.beta_radians,
    lambda = source.mass_sheet_lambda;
  if (!(zl > 0 && zs > zl && zs <= 20 && sigma > 0 &&
        sigma <= speed_of_light_m_per_s / 1000 * .01 && beta > 0 &&
        beta <= .1 && source.source_flux > 0 && source.psf_sigma_radians >= 0 &&
        source.psf_sigma_radians <= .1 && lambda >= 1e-6 && lambda <= 1)) {
    out.status_ = S::outside_domain;
    return out;
  }
  double redshifts[]{source.lens_redshift, source.source_redshift};
  p.geometry.maximum_native_bytes =
      std::min(p.geometry.maximum_native_bytes, p.maximum_native_bytes);
  cosmology::EarlyLateBatch batch;
  try {
    batch = cosmology::evaluate_early_late(
        source.background, redshifts,
        cosmology::early_late_mask(cosmology::EarlyLateOutput::dm_mpc),
        p.geometry);
  } catch (const std::bad_alloc &) {
    out.status_ = S::work_limit;
    return out;
  }
  out.callbacks_ = batch.callbacks;
  if (batch.status != S::ok || batch.rows.size() != 2) {
    out.status_ = batch.status;
    return out;
  }
  const auto &l =
      batch.rows[0].outputs[(unsigned)cosmology::EarlyLateOutput::dm_mpc];
  const auto &s =
      batch.rows[1].outputs[(unsigned)cosmology::EarlyLateOutput::dm_mpc];
  if (!l.value || !s.value) {
    out.status_ = !l.value ? l.status : s.status;
    return out;
  }
  W cl = *l.value, cs = *s.value, el = l.error_estimate, es = s.error_estimate,
    gap = cs - cl, egap = es + el + 64 * eps * (cs + cl);
  if (!(cl > el && cs > es && gap > egap)) {
    out.status_ = S::conditioning_budget_exceeded;
    return out;
  }
  W dd = cl / (1 + zl), ds = cs / (1 + zs), dds = gap / (1 + zs),
    ratio = gap / cs, ratio_low = (gap - egap) / (cs + es),
    ratio_high = (gap + egap) / (cs - es),
    ratio_error = std::max(ratio - ratio_low, ratio_high - ratio),
    amplitude = 4 * std::numbers::pi_v<W> *
                std::pow(sigma / (speed_of_light_m_per_s / 1000), 2),
    e = amplitude * ratio, ee = amplitude * ratio_error + 64 * eps * e;
  // Topology/near-caustic admission must resolve both branches at the inherited
  // distance quality, even when the caller only asks for a point-source flux.
  if (beta / lambda + e + ee > .1L) {
    out.status_ = S::outside_domain;
    return out;
  }
  W scaled = lambda * e;
  if (beta < 1e-6L * scaled + lambda * ee ||
      std::abs(beta - scaled) < 1e-6L * scaled + lambda * ee) {
    out.status_ = S::conditioning_budget_exceeded;
    return out;
  }
  project(out.geometry_[0], dd, el / (1 + zl), 1e-8);
  project(out.geometry_[1], ds, es / (1 + zs), 1e-8);
  project(out.geometry_[2], dds, egap / (1 + zs), 1e-8);
  project(out.geometry_[3], e, ee, 1e-8);
  for (const auto &v : out.geometry_)
    if (!v.value) {
      out.status_ = v.status;
      return out;
    }
  out.count_ = beta < scaled ? 2 : 1;
  for (size_t i = 0; i < out.count_; ++i) {
    W sign = i ? -1 : 1;
    out.x_[i] = beta / lambda + sign * e;
    out.x_error_[i] = ee;
    out.mu_[i] = (1 + sign * scaled / beta) / (lambda * lambda);
    out.mu_error_[i] = ee / (lambda * beta) + 64 * eps * std::abs(out.mu_[i]);
    out.flux_[i] = source.source_flux * std::abs(out.mu_[i]);
    out.flux_error_[i] = source.source_flux * out.mu_error_[i];
  }
  W coefficient = cl * cs / gap, low = (cl - el) * (cs - es) / (gap + egap),
    high = (cl + el) * (cs + es) / (gap - egap),
    coefficient_error = std::max(coefficient - low, high - coefficient),
    unit = megaparsec_in_metres_wide() / speed_of_light_m_per_s;
  out.delay_ = coefficient * 2 * e * beta * unit;
  out.delay_error_ =
      std::max(out.delay_ - (coefficient - coefficient_error) * 2 * (e - ee) *
                                beta * unit,
               (coefficient + coefficient_error) * 2 * (e + ee) * beta * unit -
                   out.delay_) +
      64 * eps * out.delay_;
  try {
    out.source_ = source;
  } catch (const std::bad_alloc &) {
    out.status_ = S::work_limit;
    return out;
  }
  out.status_ = S::ok;
  return out;
}
SISPrediction SISThinLens::predict(unsigned outputs, double r) const {
  SISPrediction out;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  if (!arithmetic() || !outputs || (outputs & ~15u) || !tolerance(r))
    return out;
  out.status = S::ok;
  out.requested = outputs;
  out.image_count = count_;
  for (size_t i = 0; i < count_; ++i) {
    auto &image = out.images[i];
    if (outputs & lens_positions)
      project(image.x_radians, x_[i], x_error_[i], r);
    if (outputs & lens_magnifications) {
      project(image.signed_magnification, mu_[i], mu_error_[i], r);
      if (image.signed_magnification.value)
        image.parity = i ? -1 : 1;
    }
    if (outputs & lens_fluxes)
      project(image.flux, flux_[i], flux_error_[i], r);
    if (outputs & lens_delays)
      project(image.relative_delay_seconds, i ? delay_ : 0,
              i ? delay_error_ : 0, r);
  }
  return out;
}
SISPixelBatch SISThinLens::pixels(std::span<const PixelRectangle> pixels,
                                  SISPixelPolicy p) const {
  SISPixelBatch out;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  if (!arithmetic() || !tolerance(p.relative_tolerance))
    return out;
  auto bytes = sis_payload_bound(pixels.size(),
                                 source_->background.drag_origin.size() + 1);
  if (pixels.size() > 65536 || pixels.size() > p.maximum_pixels || !bytes ||
      *bytes > p.maximum_native_bytes || *bytes > size_t(1024) * 1024 * 1024) {
    out.status = S::work_limit;
    return out;
  }
  try {
    out.rows.resize(pixels.size());
  } catch (const std::bad_alloc &) {
    out.status = S::work_limit;
    return out;
  }
  out.status = S::ok;
  W sigma = source_->psf_sigma_radians;
  for (size_t n = 0; n < pixels.size(); ++n) {
    const auto &pixel = pixels[n];
    auto &row = out.rows[n];
    bool finite = std::isfinite(pixel.x_low) && std::isfinite(pixel.x_high) &&
                  std::isfinite(pixel.y_low) && std::isfinite(pixel.y_high);
    if (!finite) {
      row.flux.status = S::nonfinite_input;
      continue;
    }
    if (!(pixel.x_low < pixel.x_high && pixel.y_low < pixel.y_high) ||
        std::max({std::abs(pixel.x_low), std::abs(pixel.x_high),
                  std::abs(pixel.y_low), std::abs(pixel.y_high)}) > .2) {
      row.flux.status = S::outside_domain;
      continue;
    }
    W flux = 0, error = 0;
    S cause = S::ok;
    for (size_t i = 0; i < count_; ++i) {
      if (out.image_pixel_terms == p.maximum_image_pixel_terms) {
        cause = S::work_limit;
        break;
      }
      ++out.image_pixel_terms;
      ++row.image_pixel_terms;
      W fraction = 0, fraction_error = 0;
      if (sigma == 0) {
        if ((std::abs(x_[i] - pixel.x_low) <= x_error_[i]) ||
            (std::abs(x_[i] - pixel.x_high) <= x_error_[i])) {
          cause = S::conditioning_budget_exceeded;
          break;
        }
        fraction = (x_[i] >= pixel.x_low && x_[i] < pixel.x_high &&
                    0 >= pixel.y_low && 0 < pixel.y_high)
                       ? 1
                       : 0;
      } else {
        W ex = 0, ey = 0,
          px = normal_interval((pixel.x_low - x_[i]) / sigma,
                               (pixel.x_high - x_[i]) / sigma, ex),
          py = normal_interval(pixel.y_low / sigma, pixel.y_high / sigma, ey);
        row.cdf_evaluations += 4;
        out.cdf_evaluations += 4;
        if (!(px > 0 && py > 0)) {
          cause = S::outside_domain;
          break;
        }
        fraction = px * py;
        // Bound each endpoint density over the retained position interval.
        // This tightens the frozen global derivative bound in distant tails.
        W dlo = std::max(0.L, std::abs(W(pixel.x_low) - x_[i]) - x_error_[i]) /
                sigma,
          dhi = std::max(0.L, std::abs(W(pixel.x_high) - x_[i]) - x_error_[i]) /
                sigma,
          location_error =
              x_error_[i] / (std::sqrt(2 * std::numbers::pi_v<W>) * sigma) *
              (std::exp(-dlo * dlo / 2) + std::exp(-dhi * dhi / 2));
        fraction_error = ex * py + ey * px + ex * ey + location_error * py;
        // Positive Gaussian marginals do not guarantee a representable joint
        // probability or diagnostic. Refuse wide-product underflow before it
        // can masquerade as the exact-zero point limit or a dropped image.
        if (!(fraction > 0) || !(fraction_error > 0)) {
          cause = S::conditioning_budget_exceeded;
          break;
        }
      }
      const W term = flux_[i] * fraction,
              term_error = flux_error_[i] * fraction +
                           (flux_[i] + flux_error_[i]) * fraction_error;
      if (fraction > 0 && (!(term > 0) || !(term_error > 0))) {
        cause = S::conditioning_budget_exceeded;
        break;
      }
      flux += term;
      error += term_error;
    }
    if (cause != S::ok)
      row.flux.status = cause;
    else
      project(row.flux, flux, error, p.relative_tolerance);
  }
  return out;
}
} // namespace irred::lensing
