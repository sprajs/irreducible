// Authored independent observer-time/frequency Gauss-Legendre reference.
// No production interpolation, radiometric or integration functions in
// reference. Exact SI definitions and system pi are shared ancestry, not
// independent physics.
#include "irred/temporal_photometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>
namespace {
namespace p = irred::photometry;
using W = long double;
using S = irred::numerics::Status;
constexpr W h = 6.62607015e-34L, c = 299792458.L;
unsigned checks = 0;
void check(bool b, const char *what) {
  ++checks;
  if (!b)
    throw std::runtime_error(what);
}
void near(const p::Outcome &v, W expected, const char *what) {
  check(v.availability == p::Availability::available && v.value.has_value(),
        what);
  if (expected > 0)
    check(*v.value > 0, what);
  else
    check(*v.value == 0, what);
  if (std::abs(W(*v.value) - expected) > 2e-12L * std::abs(expected) + 1e-300L)
    std::cerr << what << " actual=" << double(*v.value)
              << " expected=" << double(expected) << '\n';
  check(std::abs(W(*v.value) - expected) <=
            2e-12L * std::abs(expected) + 1e-300L,
        what);
}
W interpolate(std::span<const double> x, std::span<const double> y, W at) {
  if (at < x.front() || at > x.back())
    return 0;
  for (std::size_t j = 1; j < x.size(); ++j)
    if (at <= x[j])
      return W(y[j - 1]) +
             (at - x[j - 1]) * (W(y[j]) - y[j - 1]) / (W(x[j]) - x[j - 1]);
  return y.back();
}
W luminosity(const p::TemporalGrid &g, W time, W wavelength) {
  if (time < g.rest_time_second.front() || time > g.rest_time_second.back() ||
      wavelength < g.rest_wavelength_metre.front() ||
      wavelength > g.rest_wavelength_metre.back())
    return 0;
  for (std::size_t i = 1; i < g.rest_time_second.size(); ++i)
    if (time <= g.rest_time_second[i]) {
      const auto n = g.rest_wavelength_metre.size();
      const W a = interpolate(
                  g.rest_wavelength_metre,
                  g.luminosity_watt_per_metre.subspan((i - 1) * n, n),
                  wavelength),
              b = interpolate(g.rest_wavelength_metre,
                              g.luminosity_watt_per_metre.subspan(i * n, n),
                              wavelength);
      return a + (time - g.rest_time_second[i - 1]) * (b - a) /
                     (W(g.rest_time_second[i]) - g.rest_time_second[i - 1]);
    }
  return 0;
}
std::vector<std::pair<W, W>> legendre(unsigned n) {
  std::vector<std::pair<W, W>> rule(n);
  for (unsigned i = 0; i < (n + 1) / 2; ++i) {
    W x = std::cos(std::numbers::pi_v<W> * (W(i) + .75L) / (W(n) + .5L));
    for (unsigned iter = 0; iter < 64; ++iter) {
      W p0 = 1, p1 = x;
      for (unsigned k = 2; k <= n; ++k) {
        W next = ((2 * W(k) - 1) * x * p1 - (W(k) - 1) * p0) / k;
        p0 = p1;
        p1 = next;
      }
      const W derivative = n * (x * p1 - p0) / (x * x - 1),
              next = x - p1 / derivative;
      if (std::abs(next - x) <= 4 * std::numeric_limits<W>::epsilon()) {
        x = next;
        break;
      }
      x = next;
    }
    W p0 = 1, p1 = x;
    for (unsigned k = 2; k <= n; ++k) {
      W next = ((2 * W(k) - 1) * x * p1 - (W(k) - 1) * p0) / k;
      p0 = p1;
      p1 = next;
    }
    const W derivative = n * (x * p1 - p0) / (x * x - 1),
            weight = 2 / ((1 - x * x) * derivative * derivative);
    rule[i] = {-x, weight};
    rule[n - 1 - i] = {x, weight};
  }
  return rule;
}
std::array<W, 3> oracle(const p::TemporalGrid &g, const p::TemporalBand &band,
                        const p::TemporalExposure &e, unsigned n) {
  const auto rule = legendre(n);
  const W r = 1 + W(e.redshift), epoch = e.source_epoch_observer_second,
          d = e.luminosity_distance_metre;
  const W lo = std::max(W(e.observer_lower_second),
                        epoch + r * g.rest_time_second.front()),
          hi = std::min(W(e.observer_upper_second),
                        epoch + r * g.rest_time_second.back());
  if (hi <= lo)
    return {};
  const W wave_lo = std::max(r * g.rest_wavelength_metre.front(),
                             W(band.passband.wavelength_metre.front())),
          wave_hi = std::min(r * g.rest_wavelength_metre.back(),
                             W(band.passband.wavelength_metre.back()));
  if (wave_hi <= wave_lo)
    return {};
  std::vector<W> time_knots{lo, hi}, wave_knots{wave_lo, wave_hi};
  for (double t : g.rest_time_second)
    if (epoch + r * t > lo && epoch + r * t < hi)
      time_knots.push_back(epoch + r * t);
  for (double x : g.rest_wavelength_metre)
    if (r * x > wave_lo && r * x < wave_hi)
      wave_knots.push_back(r * x);
  for (double x : band.passband.wavelength_metre)
    if (x > wave_lo && x < wave_hi)
      wave_knots.push_back(x);
  for (auto *v : {&time_knots, &wave_knots}) {
    std::sort(v->begin(), v->end());
    v->erase(std::unique(v->begin(), v->end()), v->end());
  }
  std::array<W, 3> result{};
  for (std::size_t i = 1; i < time_knots.size(); ++i) {
    const W half_time = (time_knots[i] - time_knots[i - 1]) / 2,
            mid_time = (time_knots[i] + time_knots[i - 1]) / 2;
    for (std::size_t j = 1; j < wave_knots.size(); ++j) {
      const W nu_lo = c / wave_knots[j], nu_hi = c / wave_knots[j - 1],
              half_nu = (nu_hi - nu_lo) / 2, mid_nu = (nu_hi + nu_lo) / 2;
      for (const auto &[xt, wt] : rule) {
        const W observer_time = mid_time + half_time * xt,
                rest_time = (observer_time - epoch) / r;
        for (const auto &[xn, wn] : rule) {
          const W nu = mid_nu + half_nu * xn, lr = c / (r * nu),
                  observed = c / nu;
          const W lnu = lr * lr * luminosity(g, rest_time, lr) / c,
                  fnu = r * lnu / (4 * std::numbers::pi_v<W> * d * d),
                  transmission =
                      interpolate(band.passband.wavelength_metre,
                                  band.passband.optical_transmission, observed),
                  weight = half_time * half_nu * wt * wn;
          result[0] += weight * fnu;
          result[1] += weight * fnu * transmission;
          result[2] += weight * fnu * transmission / (h * nu);
        }
      }
    }
  }
  result[0] /= W(e.observer_upper_second) - e.observer_lower_second;
  result[1] *= e.collecting_area_square_metre;
  result[2] *= e.collecting_area_square_metre;
  return result;
}
void compare(const p::TemporalGrid &g, const p::TemporalBand &b,
             const p::TemporalExposure &e) {
  auto coarse = oracle(g, b, e, 32), fine = oracle(g, b, e, 64);
  for (unsigned i = 0; i < 3; ++i) {
    if (fine[i] > 0)
      check(std::abs(coarse[i] - fine[i]) <= 2e-13L * std::abs(fine[i]),
            "independent time/frequency refinement allocation");
    else
      check(coarse[i] == 0, "zero reference stable");
  }
  std::array grids{g};
  std::array bands{b};
  auto prepared = p::prepare_temporal(grids, bands);
  check(prepared.status() == S::ok, "peer preparation");
  std::array exposures{e};
  auto out = p::evaluate_temporal(prepared, exposures);
  check(out.status == S::ok && out.rows[0].admission_status == S::ok,
        "peer exposure admission");
  near(out.rows[0].mean_flux_watt_per_square_metre, fine[0],
       "independent frequency-time mean flux");
  near(out.rows[0].energy_joule, fine[1], "independent frequency-time energy");
  near(out.rows[0].expected_photons, fine[2],
       "independent frequency-time photons");
}
} // namespace
int main() {
  try {
    constexpr double m = 0x1p-20;
    std::array t{0., 1., 2.};
    std::array x{m, 3 * m};
    std::array l{7., 17., 17., 41., 27., 65.};
    std::array bw{2 * m, 4 * m};
    std::array trans{.25, .75};
    // L(τ,X)=2+3τ+5X+7τX, X=λ_rest/m; original bilinear polynomial.
    p::TemporalGrid grid{"bilinear-polynomial",
                         "independent authored bilinear synthetic control", t,
                         x, l};
    p::TemporalBand band{"linear-band",
                         "declared fixed synthetic optical transmission",
                         {bw, trans}};
    p::TemporalExposure e{0, 0, 2, 1, 2, 0, -2, 6};
    std::array grids{grid};
    std::array bands{band};
    auto prepared = p::prepare_temporal(grids, bands);
    std::array exposures{e};
    auto out = p::evaluate_temporal(prepared, exposures);
    check(out.rows[0].admission_status == S::ok &&
              out.rows[0].coverage == p::TemporalCoverage::partial &&
              out.rows[0].covered_fraction == .5L,
          "bilinear partial support");
    near(out.rows[0].mean_flux_watt_per_square_metre,
         23.L * m / (32 * std::numbers::pi_v<W>),
         "bilinear original polynomial incident integral");
    near(out.rows[0].energy_joule, 6.L * m / std::numbers::pi_v<W>,
         "bilinear original polynomial collected integral");
    near(out.rows[0].expected_photons,
         467.L / 24 * m * m / (std::numbers::pi_v<W> * h * c),
         "bilinear original polynomial photon integral");
    compare(grid, band, e);
    e.observer_lower_second = .5;
    e.observer_upper_second = 3.5;
    compare(grid, band, e);
    e.source_epoch_observer_second = 100;
    e.observer_lower_second = 100.5;
    e.observer_upper_second = 103.5;
    compare(grid, band, e);
    // Exactly linear temporal and spectral knot refinement changes no model.
    std::array rt{0., .5, 1., 1.5, 2.};
    std::array rx{m, 2 * m, 3 * m};
    std::array<double, 15> rl{};
    for (std::size_t i = 0; i < rt.size(); ++i)
      for (std::size_t j = 0; j < rx.size(); ++j) {
        const W X = W(rx[j]) / m;
        rl[i * rx.size() + j] = 2 + 3 * rt[i] + 5 * X + 7 * rt[i] * X;
      }
    p::TemporalGrid refined{"refined-grid", "exactly linear inserted knots", rt,
                            rx, rl};
    e = {0, 0, 2, 1, 2, 0, -2, 6};
    std::array rgrids{refined};
    auto rprepared = p::prepare_temporal(rgrids, bands);
    auto refined_out = p::evaluate_temporal(rprepared, std::span(&e, 1));
    near(refined_out.rows[0].mean_flux_watt_per_square_metre,
         *out.rows[0].mean_flux_watt_per_square_metre.value,
         "temporal/spectral knot refinement mean");
    near(refined_out.rows[0].energy_joule, *out.rows[0].energy_joule.value,
         "temporal/spectral knot refinement energy");
    near(refined_out.rows[0].expected_photons,
         *out.rows[0].expected_photons.value,
         "temporal/spectral knot refinement photons");
    compare(refined, band, e);
    std::array time2{-2., 0., 1., 4.};
    std::array wave2{m, 1.5 * m, 2.25 * m, 4 * m};
    std::array lum2{1., 3., 2., 4., 2., 8., .5, 3.,
                    4., 7., 1., 6., 2., .2, 8., 1.};
    std::array bandwave{1.2 * m, 1.75 * m, 2 * m, 2.8 * m, 3.7 * m};
    std::array bandt{0., .8, .1, 1., .3};
    p::TemporalGrid mixed{"mixed-grid",
                          "independent nonseparable synthetic grid", time2,
                          wave2, lum2};
    p::TemporalBand mixed_band{
        "mixed-band", "fixed irregular synthetic response", {bandwave, bandt}};
    e = {0, 0, 1e12, .7, 2, 8, 6, 11};
    compare(mixed, mixed_band, e);
    e.observer_lower_second = 0;
    e.observer_upper_second = 20;
    compare(mixed, mixed_band, e);
    e.redshift = 0;
    e.observer_lower_second = 7;
    e.observer_upper_second = 10;
    compare(mixed, mixed_band, e);
    // Domain admission precedes no-overlap/zero shortcuts.
    bandt.fill(0);
    compare(mixed, mixed_band, e);
    lum2[0] = -1;
    std::array mg{mixed};
    std::array mb{mixed_band};
    check(p::prepare_temporal(mg, mb).status() == S::outside_domain,
          "zero transmission cannot hide invalid source");
    std::cout << "PASS " << checks
              << " temporal photometry independent peer controls GL32/64\n";
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
