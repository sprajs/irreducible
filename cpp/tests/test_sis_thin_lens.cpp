// Original synthetic controls: elementary radiation-only chi, stationarity of
// the explicit potential, radial/tangential Jacobian, Fermat difference, and
// fixed-panel GL8 integration of a Gaussian density (not production erf/CDF).
// Angles/mu/flux allocate2e-12 relative; delay2e-11; pixels2e-12 absolute
// source-flux units. 64/128 panels must refine within5% of the pixel
// allocation.
#include "irred/sis_thin_lens.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred::lensing;
using S = irred::numerics::Status;
using W = long double;
namespace {
unsigned checks = 0;
W maximum_pixel_difference = 0, maximum_reference_fraction = 0;
void need(bool b, const char *label) {
  ++checks;
  if (!b)
    throw std::runtime_error(label);
}
void near(W a, W b, W absolute, W relative, const char *label) {
  if (std::abs(a - b) > absolute + relative * std::abs(b))
    std::cerr << label << " got=" << (double)a << " ref=" << (double)b
              << " difference=" << (double)std::abs(a - b) << '\n';
  need(std::abs(a - b) <= absolute + relative * std::abs(b), label);
}
W value(const LensValue &v) {
  need(v.status == S::ok && v.value.has_value(), "accepted scalar");
  return *v.value;
}
W theta_e() {
  return 2 * std::numbers::pi_v<W> * std::pow(220 / 299792.458L, 2);
}
SISSource source() {
  return {{{70, 0, 1, 0, 1}, 1059, "synthetic radiation-only geometry"},
          .5,
          2,
          220,
          (double)(.4L * theta_e()),
          1,
          (double)(.2L * theta_e()),
          1};
}
W chi(W z, W h = 70) { return 299792.458L / h * z / (1 + z); }
W potential(W x, W e, W lambda) {
  return lambda * e * std::abs(x) + (1 - lambda) * x * x / 2;
}
W phi(W x, W beta, W e, W lambda) {
  return (x - beta) * (x - beta) / 2 - potential(x, e, lambda);
}
W gaussian(W lo, W hi, unsigned panels) {
  const W x[]{
      .183434642495649804939476142360184L, .525532409916328985817739049189246L,
      .796666477413626739591553936475831L, .960289856497536231683560868569473L};
  const W weight[]{
      .362683783378361982965150449277195L, .313706645877887287337962201986601L,
      .222381034453374470544355994426240L, .101228536290376259152531354309962L};
  W sum = 0, h = (hi - lo) / panels;
  for (unsigned i = 0; i < panels; ++i) {
    W middle = lo + (i + .5L) * h;
    for (unsigned k = 0; k < 4; ++k)
      for (W sign : {-1.L, 1.L}) {
        W t = middle + sign * x[k] * h / 2;
        sum += weight[k] * std::exp(-t * t / 2) * h / 2 /
               std::sqrt(2 * std::numbers::pi_v<W>);
      }
  }
  return sum;
}
W pixel_reference(const SISSource &s, const PixelRectangle &p,
                  unsigned panels) {
  W e = theta_e(), beta = s.beta_radians, lambda = s.mass_sheet_lambda,
    sigma = s.psf_sigma_radians, total = 0;
  for (W sign : {1.L, -1.L}) {
    W x = beta / lambda + sign * e;
    W mu = (1 + sign * lambda * e / beta) / (lambda * lambda);
    total += s.source_flux * std::abs(mu) *
             gaussian((p.x_low - x) / sigma, (p.x_high - x) / sigma, panels) *
             gaussian(p.y_low / sigma, p.y_high / sigma, panels);
  }
  return total;
}
void image_controls(const SISSource &s, const SISThinLens &lens) {
  auto prediction = lens.predict(15);
  need(prediction.image_count == 2, "two images");
  W e = theta_e(), beta = s.beta_radians, lambda = s.mass_sheet_lambda;
  W cl = chi(s.lens_redshift, s.background.model.h0_km_s_mpc),
    cs = chi(s.source_redshift, s.background.model.h0_km_s_mpc);
  near(value(lens.geometry()[0]), cl / (1 + s.lens_redshift), 1e-9, 2e-11,
       "Dd analytic");
  near(value(lens.geometry()[1]), cs / (1 + s.source_redshift), 1e-9, 2e-11,
       "Ds analytic");
  near(value(lens.geometry()[2]), (cs - cl) / (1 + s.source_redshift), 1e-9,
       2e-11, "Dds analytic");
  near(value(lens.geometry()[3]), e, 0, 2e-12, "Einstein angle analytic");
  W x[2]{};
  for (unsigned i = 0; i < 2; ++i) {
    auto &image = prediction.images[i];
    W sign = i ? -1 : 1;
    x[i] = value(image.x_radians);
    near(x[i], beta / lambda + sign * e, 0, 2e-12, "image analytic");
    // Differentiate lens map beta=lambda*(theta-e theta/|theta|): radial
    // eigenvalue lambda, tangential lambda*(1-e/|theta|), mu=1/determinant.
    W jacobian_mu = 1 / (lambda * lambda * (1 - e / std::abs(x[i])));
    near(value(image.signed_magnification), jacobian_mu, 0, 2e-12,
         "Jacobian magnification");
    need(image.parity == (i ? -1 : 1), "signed image parity");
    near(value(image.flux), s.source_flux * std::abs(jacobian_mu), 0, 2e-12,
         "image flux analytic");
    W stationarity = x[i] - beta - (lambda * e * sign + (1 - lambda) * x[i]);
    near(stationarity, 0, 2e-12 * e, 0, "Fermat stationary image");
  }
  // Original IAU exact pc definition; no production constant/helper readback.
  W mpc = 1e6L * (648000 / std::numbers::pi_v<W>)*149597870700.L;
  W delay = (phi(x[1], beta, e, lambda) - phi(x[0], beta, e, lambda)) * cl *
            cs / (cs - cl) * mpc / 299792458;
  near(value(prediction.images[0].relative_delay_seconds), 0, 0, 0,
       "earliest arrival zero");
  near(value(prediction.images[1].relative_delay_seconds), delay, 0, 2e-11,
       "independent Fermat delay");
}
} // namespace
int main() {
  try {
    auto s = source();
    auto lens = prepare_sis_thin_lens(s);
    need(lens.status() == S::ok, "prepared lens");
    need(lens.preparation_callbacks() > 0, "geometry callback accounting");
    image_controls(s, lens);
    W e = theta_e();
    std::vector<PixelRectangle> pixels;
    for (int iy = -5; iy < 5; ++iy)
      for (int ix = -12; ix < 12; ++ix)
        pixels.push_back({(double)(ix * .2L * e), (double)((ix + 1) * .2L * e),
                          (double)(iy * .2L * e),
                          (double)((iy + 1) * .2L * e)});
    auto result = lens.pixels(pixels);
    need(result.rows.size() == pixels.size(), "ordered pixels");
    need(result.image_pixel_terms == 2 * pixels.size() &&
             result.cdf_evaluations == 8 * pixels.size(),
         "actual pixel and CDF counts");
    for (size_t i = 0; i < pixels.size(); ++i) {
      W ref = pixel_reference(s, pixels[i], 128),
        coarse = pixel_reference(s, pixels[i], 64);
      W fraction = std::abs(ref - coarse) / 2e-12L;
      maximum_reference_fraction =
          std::max(maximum_reference_fraction, fraction);
      near(ref, coarse, .05L * 2e-12, 0, "Gaussian refinement");
      W got = value(result.rows[i].flux);
      maximum_pixel_difference =
          std::max(maximum_pixel_difference, std::abs(got - ref));
      near(got, ref, 2e-12, 0, "independent Gaussian pixel");
    }
    // Larger finite partition: omitted tails beyond8sigma have negligible flux
    // relative to2e-12. Compare full photon/linear-flux sum with sum |mu|
    // Fsource.
    std::vector<PixelRectangle> partition;
    for (int iy = -8; iy < 8; ++iy)
      for (int ix = -15; ix < 15; ++ix)
        partition.push_back(
            {(double)(ix * .2L * e), (double)((ix + 1) * .2L * e),
             (double)(iy * .2L * e), (double)((iy + 1) * .2L * e)});
    auto integrated = lens.pixels(partition);
    W sum = 0;
    for (auto &row : integrated.rows)
      sum += value(row.flux);
    auto prediction = lens.predict(15);
    near(sum,
         value(prediction.images[0].flux) + value(prediction.images[1].flux),
         2e-12, 0, "finite field flux conservation");
    // Synthetic recovery conditional on SIS lambda=1: positions fix E and beta,
    // image-flux difference fixes source flux, and delay fixes distance
    // coefficient.
    W xp = value(prediction.images[0].x_radians),
      xm = value(prediction.images[1].x_radians),
      recovered_beta = (xp + xm) / 2, recovered_e = (xp - xm) / 2,
      recovered_flux = (value(prediction.images[0].flux) -
                        value(prediction.images[1].flux)) /
                       2,
      recovered_coefficient =
          value(prediction.images[1].relative_delay_seconds) /
          (2 * recovered_e * recovered_beta);
    near(recovered_beta, s.beta_radians, 0, 2e-12, "synthetic source recovery");
    near(recovered_e, e, 0, 2e-12, "synthetic template recovery");
    near(recovered_flux, s.source_flux, 0, 2e-12,
         "synthetic source flux recovery");
    W mpc = 1e6L * 648000 / std::numbers::pi_v<W> * 149597870700.L;
    near(recovered_coefficient,
         chi(.5) * chi(2) / (chi(2) - chi(.5)) * mpc / 299792458, 0, 2e-11,
         "synthetic delay recovery");
    for (double lambda : {.2, .5, .9}) {
      auto transformed = s;
      transformed.mass_sheet_lambda = lambda;
      transformed.beta_radians *= lambda;
      transformed.source_flux *= lambda * lambda;
      auto family = prepare_sis_thin_lens(transformed);
      image_controls(transformed, family);
      auto p = family.predict(15);
      auto f = family.pixels(pixels);
      for (unsigned i = 0; i < 2; ++i) {
        near(value(p.images[i].x_radians),
             value(prediction.images[i].x_radians), 0, 2e-12,
             "MST same images");
        near(value(p.images[i].flux), value(prediction.images[i].flux), 0,
             2e-12, "MST same flux");
      }
      near(value(p.images[1].relative_delay_seconds),
           lambda * value(prediction.images[1].relative_delay_seconds), 0,
           2e-11, "MST scaled delay");
      for (size_t i = 0; i < pixels.size(); ++i)
        near(value(f.rows[i].flux), value(result.rows[i].flux), 2e-12, 0,
             "MST same pixels");
      transformed.background.model.h0_km_s_mpc *= lambda;
      auto degenerate = prepare_sis_thin_lens(transformed).predict(15);
      near(value(degenerate.images[1].relative_delay_seconds),
           value(prediction.images[1].relative_delay_seconds), 0, 2e-11,
           "MST H0 joint delay degeneracy");
    }
    auto scaled = s;
    scaled.background.model.h0_km_s_mpc *= 2;
    auto h0 = prepare_sis_thin_lens(scaled).predict(15);
    near(value(h0.images[0].x_radians), xp, 0, 2e-12,
         "fixed fraction H0 image scaling");
    near(value(h0.images[1].relative_delay_seconds),
         value(prediction.images[1].relative_delay_seconds) / 2, 0, 2e-11,
         "fixed fraction H0 delay scaling");
    auto single = s;
    single.beta_radians = (double)(2 * e);
    need(prepare_sis_thin_lens(single).predict(15).image_count == 1,
         "one image outside SIS ring");
    auto zero = s;
    zero.psf_sigma_radians = 0;
    auto points = prepare_sis_thin_lens(zero);
    PixelRectangle box[]{{(double)(-2 * e), (double)(2 * e), 0, (double)e},
                         {(double)(-2 * e), (double)(2 * e), (double)-e, 0}};
    auto point_result = points.pixels(box);
    near(value(point_result.rows[0].flux), 5, 0, 2e-12,
         "zero PSF halfopen y boundary");
    near(value(point_result.rows[1].flux), 0, 0, 0,
         "zero PSF complementary partition");
    need(point_result.cdf_evaluations == 0 &&
             point_result.image_pixel_terms == 4,
         "point limit counts");
    PixelRectangle mixed[]{pixels[0],
                           {1, 0, 0, 1},
                           {0, 1, 0, std::numeric_limits<double>::quiet_NaN()},
                           pixels[1]};
    auto mixed_result = lens.pixels(mixed);
    need(mixed_result.rows[0].flux.value && !mixed_result.rows[1].flux.value &&
             !mixed_result.rows[2].flux.value &&
             mixed_result.rows[3].flux.value,
         "mixed pixel states");
    SISPixelPolicy quota;
    quota.maximum_image_pixel_terms = 2;
    auto capped = lens.pixels(std::span(pixels.data(), 2), quota);
    need(capped.rows[0].flux.value && !capped.rows[1].flux.value &&
             capped.image_pixel_terms == 2,
         "pixel work cap partial states");
    quota = {};
    quota.maximum_native_bytes = 0;
    need(lens.pixels(pixels, quota).rows.empty(), "pixel payload preflight");
    quota = {};
    quota.maximum_pixels = 0;
    need(lens.pixels(pixels, quota).rows.empty(), "pixel count preflight");
    quota = {};
    quota.relative_tolerance = 1e-30;
    need(!lens.pixels(std::span(pixels.data(), 1), quota).rows[0].flux.value &&
             lens.predict(1).images[0].x_radians.value,
         "pixel refusal preserves image prediction");
    need(!sis_payload_bound(SIZE_MAX, 1), "payload arithmetic overflow");
    for (unsigned mask : {1u, 2u, 4u, 8u}) {
      auto p = lens.predict(mask);
      need(bool(p.images[0].x_radians.value) == bool(mask & 1) &&
               bool(p.images[0].signed_magnification.value) == bool(mask & 2) &&
               bool(p.images[0].flux.value) == bool(mask & 4) &&
               bool(p.images[0].relative_delay_seconds.value) == bool(mask & 8),
           "requested groups only");
    }
    SISPreparationPolicy prep;
    prep.geometry.maximum_total_callbacks = 0;
    need(prepare_sis_thin_lens(s, prep).status() == S::work_limit,
         "geometry callback cap");
    prep = {};
    prep.maximum_native_bytes = 1;
    need(prepare_sis_thin_lens(s, prep).status() == S::work_limit,
         "preparation payload cap");
    for (unsigned k = 0; k < 8; ++k) {
      auto invalid = s;
      if (k == 0)
        invalid.lens_redshift = invalid.source_redshift;
      if (k == 1)
        invalid.beta_radians = 0;
      if (k == 2)
        invalid.beta_radians = (double)e;
      if (k == 3)
        invalid.mass_sheet_lambda = 1.01;
      if (k == 4)
        invalid.psf_sigma_radians = -1;
      if (k == 5)
        invalid.background.model.omega_r = 0;
      if (k == 6)
        invalid.source_redshift = std::nextafter(invalid.lens_redshift, 1.);
      if (k == 7)
        invalid.sigma_v_km_s = std::numeric_limits<double>::infinity();
      need(prepare_sis_thin_lens(invalid).status() != S::ok,
           "invalid domain or unresolved geometry");
    }
    // Projection failures remain separate: extreme source flux must not erase
    // finite geometry, magnification, or arrival predictions.
    auto bright = source();
    bright.source_flux = std::numeric_limits<double>::max();
    auto bright_prediction = prepare_sis_thin_lens(bright).predict(15);
    need(bright_prediction.images[0].x_radians.value &&
             bright_prediction.images[0].signed_magnification.value &&
             !bright_prediction.images[0].flux.value &&
             bright_prediction.images[1].relative_delay_seconds.value,
         "flux overflow preserves other image groups");
    PixelRectangle boundary{(double)xp, (double)(2 * e), (double)-e, (double)e};
    need(!points.pixels(std::span(&boundary, 1)).rows[0].flux.value,
         "zero PSF unresolved x boundary refuses");
    PixelRectangle distant{.01, .02, .01, .02};
    need(!lens.pixels(std::span(&distant, 1)).rows[0].flux.value,
         "positive Gaussian underflow refuses rather than zero");
    // Both Gaussian marginals remain positive in wide arithmetic, but their
    // joint probability is below its exponent range. The former implementation
    // exposed an available flux0/error0 here. Preserve a typed refusal rather
    // than changing the unchanged analytic/GL comparison allocation.
    const double width = s.psf_sigma_radians;
    const PixelRectangle positive_product{115 * width, 116 * width,
                                          108 * width, 109 * width};
    const auto joint = lens.pixels(std::span(&positive_product, 1));
    need(joint.rows[0].flux.status == S::conditioning_budget_exceeded &&
             !joint.rows[0].flux.value && joint.image_pixel_terms == 1 &&
             joint.cdf_evaluations == 4,
         "positive Gaussian joint underflow refuses with attempted work");
    // A representable wide probability can still lose its positive weighted
    // term before summation. The same rule excludes a silently dropped image.
    auto dim = source();
    dim.source_flux = std::numeric_limits<double>::denorm_min();
    const auto dim_lens = prepare_sis_thin_lens(dim);
    const PixelRectangle positive_weighted{112 * width, 113 * width,
                                           105 * width, 106 * width};
    const auto weighted = dim_lens.pixels(std::span(&positive_weighted, 1));
    need(weighted.rows[0].flux.status == S::conditioning_budget_exceeded &&
             !weighted.rows[0].flux.value,
         "positive Gaussian weighted-term underflow refuses");
    need(lens.predict(0).image_count == 0 &&
             lens.predict(16).image_count == 0 &&
             lens.predict(1, std::numeric_limits<double>::quiet_NaN())
                     .image_count == 0,
         "invalid projection masks and budgets");
    need(lens.pixels({}).status == S::ok && lens.pixels({}).rows.empty(),
         "empty pixel batch explicit success");
    auto reversed = pixels;
    std::reverse(reversed.begin(), reversed.end());
    auto reversed_result = lens.pixels(reversed);
    for (size_t i = 0; i < pixels.size(); ++i)
      need(reversed_result.rows[i].flux.value ==
               result.rows[pixels.size() - 1 - i].flux.value,
           "pixel order permutation");
    quota = {};
    quota.maximum_image_pixel_terms = 1;
    auto incomplete = lens.pixels(std::span(pixels.data(), 1), quota);
    need(!incomplete.rows[0].flux.value && incomplete.image_pixel_terms == 1 &&
             incomplete.cdf_evaluations == 4,
         "partial image sum refuses rather than dropping image");
    auto unbounded_angle = source();
    unbounded_angle.beta_radians = .1;
    unbounded_angle.mass_sheet_lambda = 1e-6;
    need(prepare_sis_thin_lens(unbounded_angle).status() != S::ok,
         "transformed image angles stay in tangent-plane domain");
    auto temporary = prepare_sis_thin_lens(source());
    need(temporary.source() && temporary.predict(1).images[0].x_radians.value,
         "source survives temporary input destruction");
    auto copy = lens;
    auto moved = std::move(lens);
    auto *self = &moved;
    moved = std::move(*self);
    s.background.drag_origin = "mutated";
    s.beta_radians = 0;
    need(copy.source()->background.drag_origin ==
                 "synthetic radiation-only geometry" &&
             copy.predict(1).images[0].x_radians.value &&
             moved.predict(1).images[0].x_radians.value &&
             lens.predict(1).image_count == 0,
         "owned immutable source lifetime");
    int mode = std::fegetround();
    need(std::fesetround(FE_UPWARD) == 0, "set hostile rounding");
    auto rejected = moved.pixels(pixels);
    need(std::fesetround(mode) == 0, "restore rounding");
    need(rejected.rows.empty(), "arithmetic refusal");
    std::cout << "PASS " << checks
              << " SIS synthetic controls; max pixel difference="
              << (double)maximum_pixel_difference
              << " reference allocation fraction="
              << (double)maximum_reference_fraction << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
