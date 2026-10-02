// Original independent Fourier density/CDF inversion with composite Simpson
// refinement. No production Poisson sum or radiometric interpolation is used.
#include "irred/optical_detector.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
using namespace irred::photometry;
using S = irred::numerics::Status;
using W = long double;
namespace {
unsigned checks = 0;
void need(bool value, const char *why) {
  ++checks;
  if (!value) throw std::runtime_error(why);
}
void near(W a, W b, W fraction = 1) {
  need(std::abs(a - b) <= fraction * 2e-10L * (1 + std::abs(b)),
       "independent predeclared log allocation/refinement");
}
template<class F> W simpson(F f, W end, unsigned intervals) {
  const W h = end / intervals;
  W sum = f(0) + f(end);
  for (unsigned i = 1; i < intervals; ++i) sum += (i & 1 ? 4 : 2) * f(i * h);
  return sum * h / 3;
}
struct Reference { W below, density; };
Reference invert(W mean, W sd, W gain, W bias, W threshold_adu,
                 W measured_adu, unsigned intervals) {
  const W threshold = gain * (threshold_adu - bias),
          observed = gain * (measured_adu - bias),
          end = 24 / sd, pi = std::numbers::pi_v<W>;
  auto amplitude = [&](W t) { return std::exp(mean * (std::cos(t) - 1) - sd * sd * t * t / 2); };
  auto cdf = [&](W t) {
    return t == 0 ? mean - threshold
                  : amplitude(t) * std::sin(mean * std::sin(t) - threshold * t) / t;
  };
  auto density = [&](W t) { return amplitude(t) * std::cos(mean * std::sin(t) - observed * t); };
  return {.5L - simpson(cdf, end, intervals) / pi,
          gain * simpson(density, end, intervals) / pi};
}
void run() {
  std::array<double, 2> wave{1, 2}, luminosity{8e-24, 8e-24};
  std::array<double, 2> a0{.2, .4}, b0{.3, .1}, a1{.8, .6}, b1{.5, .9};
  std::array<CalibrationBand, 2> bands{{{"a", wave, 1, 2}, {"b", wave, 1, 3}}};
  std::array<std::span<const double>, 2> t0{a0, b0}, t1{a1, b1};
  std::array<CalibrationState, 2> states{{{"s0", 2, t0}, {"s1", 5, t1}}};
  std::array<DetectorBand, 2> detector{{
      {irred::detector::PhotonLaw::poisson_arrivals, .8, .3, .2, .75, 2, .5},
      {irred::detector::PhotonLaw::poisson_arrivals, .45, .7, .1, 1.5, .5, -.25}}};
  std::array<irred::detector::Observation, 2> detected{{{2, true, {}, 3.5}, {3, true, {}, 7}}};
  std::array<irred::detector::Observation, 2> censored{{{2, false, {}, {}}, {3, true, {}, 7}}};
  std::array<OpticalDetectorRecord, 3> records{{
      {"joint", detected}, {"censored", censored},
      {"selected", detected, irred::detector::SelectionMeasure::selected_only}}};
  OpticalDetectorInput input{{{wave, luminosity}, 1, 0, bands, states,
                              "synthetic-constant-spectrum", "linear-optical-states",
                              "relative-masses-2-5", "shared-state"},
                             detector, records, ConditionalDetectorLaw::independent_poisson_and_read,
                             "synthetic-detectors", "synthetic-censored-records", "conditional-independent"};
  auto result = evaluate_optical_detector(input);
  need(result.status == S::ok && result.attempts.size() == 4 && result.records.size() == 3,
       "continuous two-state two-band SDK law");
  std::array<W, 3> coarse{}, fine{};
  W coarse_event = 0, fine_event = 0, marginal_a = 0, marginal_b = 0;
  for (unsigned s = 0; s < 2; ++s) {
    const W mass = W(states[s].relative_mass) / 7;
    std::array<Reference, 2> refs_coarse{}, refs_fine{};
    for (unsigned b = 0; b < 2; ++b) {
      // Linear T through wavelengths 1,2: integral lambda*T(lambda) d lambda
      // = 2*T(1)/3 + 5*T(2)/6, independently derived polynomial antiderivative.
      const auto transmission = states[s].optical_transmission[b];
      const W integral = (2 * W(transmission[0])) / 3 + (5 * W(transmission[1])) / 6;
      const W photons = W(luminosity[0]) * bands[b].observer_exposure_second * integral /
                        (4 * std::numbers::pi_v<W> * 6.62607015e-34L * 299792458.L);
      const auto &d = detector[b];
      const W mean = photons * d.quantum_efficiency + d.background_expected_electrons +
                     W(d.dark_electrons_per_second) * bands[b].observer_exposure_second;
      refs_coarse[b] = invert(mean, d.read_noise_rms_electrons, d.gain_electrons_per_adu,
                              d.bias_adu, detected[b].detection_threshold_adu,
                              *detected[b].measured_adu, 8192);
      refs_fine[b] = invert(mean, d.read_noise_rms_electrons, d.gain_electrons_per_adu,
                            d.bias_adu, detected[b].detection_threshold_adu,
                            *detected[b].measured_adu, 16384);
      need(refs_fine[b].below > 0 && refs_fine[b].below < 1 && refs_fine[b].density > 0,
           "positive independently inverted density/CDF");
      near(std::log(refs_coarse[b].below), std::log(refs_fine[b].below), .05L);
      near(std::log(refs_coarse[b].density), std::log(refs_fine[b].density), .05L);
      near(*result.attempts[2 * s + b].nominal.rows[0].log_value,
           std::log(refs_fine[b].density));
      near(*result.attempts[2 * s + b].nominal.rows[3].log_value,
           std::log(refs_fine[b].below));
      need(std::abs(W(result.attempts[2 * s + b].nominal.source.expected_transmitted_photons) - photons)
               <= 2e-12L * photons, "quantum efficiency follows optical photons");
    }
    coarse[0] += mass * refs_coarse[0].density * refs_coarse[1].density;
    fine[0] += mass * refs_fine[0].density * refs_fine[1].density;
    coarse[1] += mass * refs_coarse[0].below * refs_coarse[1].density;
    fine[1] += mass * refs_fine[0].below * refs_fine[1].density;
    coarse_event += mass * (1 - refs_coarse[0].below) * (1 - refs_coarse[1].below);
    fine_event += mass * (1 - refs_fine[0].below) * (1 - refs_fine[1].below);
    marginal_a += mass * (1 - refs_fine[0].below);
    marginal_b += mass * (1 - refs_fine[1].below);
    near(*result.conditional[3 * s].log_record,
         std::log(refs_fine[0].density * refs_fine[1].density));
    near(*result.conditional[3 * s + 1].log_record,
         std::log(refs_fine[0].below * refs_fine[1].density));
  }
  coarse[2] = coarse[0] / coarse_event;
  fine[2] = fine[0] / fine_event;
  near(std::log(coarse_event), std::log(fine_event), .05L);
  need(std::abs(fine_event - marginal_a * marginal_b) > 1e-4L,
       "selected denominator is same joint mixture event rather than marginal product");
  for (unsigned r = 0; r < 3; ++r) {
    near(std::log(coarse[r]), std::log(fine[r]), .05L);
    near(*result.records[r].log_value, std::log(fine[r]));
    near(*result.records[r].log_all_band_detection, std::log(fine_event));
    need(result.records[r].numerical_error_estimate > 0,
         "optical and detector empirical diagnostics carried to mixture");
  }
  // A separate zero-photon Gaussian limit has no Poisson inversion ancestry.
  luminosity = {0, 0};
  for (auto &d : detector) { d.background_expected_electrons = 0; d.dark_electrons_per_second = 0; }
  auto zero = evaluate_optical_detector(input);
  need(zero.status == S::ok, "zero-photon Gaussian limit");
  W direct_density = 1, direct_event = 1;
  for (unsigned b = 0; b < 2; ++b) {
    const auto &d = detector[b];
    const W z = W(d.gain_electrons_per_adu) * (*detected[b].measured_adu - d.bias_adu) / d.read_noise_rms_electrons;
    const W t = W(d.gain_electrons_per_adu) * (detected[b].detection_threshold_adu - d.bias_adu) / d.read_noise_rms_electrons;
    direct_density *= d.gain_electrons_per_adu * std::exp(-z * z / 2) /
                      (d.read_noise_rms_electrons * std::sqrt(2 * std::numbers::pi_v<W>));
    direct_event *= std::erfc(t / std::sqrt(2.L)) / 2;
  }
  near(*zero.records[0].log_value, std::log(direct_density));
  near(*zero.records[2].log_value, std::log(direct_density / direct_event));
}
} // namespace
int main() {
  try {
    run();
    std::cout << "PASS " << checks << " independent optical-mixture Fourier controls\n";
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n'; return 1;
  }
}
