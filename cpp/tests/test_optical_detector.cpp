#include "irred/optical_detector.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
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
void near(W value, W reference) {
  need(std::isfinite(value) && std::abs(value - reference) <= 2e-10L * (1 + std::abs(reference)),
       "frozen independent log allocation");
}
W poisson(W mean, unsigned count) {
  W p = std::exp(-mean);
  for (unsigned n = 1; n <= count; ++n) p *= mean / n;
  return p;
}
struct Fixture {
  std::array<double, 2> wave{1, 2}, luminosity{8e-24, 8e-24}, low{.25, .25}, high{.75, .75};
  std::array<CalibrationBand, 2> bands{{{"blue", wave, 1, 1}, {"red", wave, 1, 1}}};
  std::array<std::span<const double>, 2> t0{low, low}, t1{high, high};
  std::array<CalibrationState, 2> states{{{"dim", 1, t0}, {"bright", 3, t1}}};
  std::array<DetectorBand, 2> detectors{{
      {irred::detector::PhotonLaw::poisson_arrivals, 1},
      {irred::detector::PhotonLaw::poisson_arrivals, .5}}};
  std::array<irred::detector::Observation, 2> observations{{{1, true, 2, {}}, {1, true, 3, {}}}};
  std::array<OpticalDetectorRecord, 1> records{{{"event", observations}}};
  OpticalDetectorInput input{{{wave, luminosity}, 1, 0, bands, states,
                              "synthetic-spectrum", "two-transmission-states",
                              "relative-masses-1-3", "same-optical-state-for-both-bands"},
                             detectors, records, ConditionalDetectorLaw::independent_poisson_and_read,
                             "fixed-synthetic-detectors", "synthetic-records", "declared-conditional-independence"};
  W full_photons() const {
    // Independent constant-spectrum wavelength polynomial, exact-binary L input.
    return W(luminosity[0]) * 1.5L /
           (4 * std::numbers::pi_v<W> * 6.62607015e-34L * 299792458.L);
  }
};
void basic() {
  Fixture f;
  auto r = evaluate_optical_detector(f.input);
  need(r.status == S::ok && r.records.size() == 1 && r.conditional.size() == 2,
       "two-state two-band SDK consumer");
  need(r.optical.aggregation == CalibrationAggregation::state_resolved_only && !r.optical.moments,
       "state resolved photometry instead of moment approximation");
  need(r.attempts.size() == 4 && r.optical.normalized_state_mass.size() == 2,
       "all required state-band attempts and masses");
  need(r.optical.normalized_state_mass[0] == .25L && r.optical.normalized_state_mass[1] == .75L,
       "single normalized optical owner");
  W joint = 0, selected = 0, marginal0 = 0, marginal1 = 0;
  for (unsigned s = 0; s < 2; ++s) {
    W n = f.full_photons() * (s ? .75L : .25L), mass = s ? .75L : .25L;
    W a = poisson(n, 2), b = poisson(n / 2, 3);
    joint += mass * a * b;
    selected += mass * (-std::expm1(-n)) * (-std::expm1(-n / 2));
    marginal0 += mass * a;
    marginal1 += mass * b;
    near(*r.conditional[s].log_record, std::log(a * b));
    near(*r.conditional[s].log_all_band_detection,
         std::log((-std::expm1(-n)) * (-std::expm1(-n / 2))));
    for (unsigned band = 0; band < 2; ++band) {
      const auto &attempt = r.attempts[2 * s + band];
      need(attempt.state_id == f.states[s].id && attempt.band_id == f.bands[band].id,
           "state then band order");
      need(attempt.nominal.rows.size() == 2 && attempt.lower_photons && attempt.upper_photons,
           "original record threshold probe and optical sensitivity endpoints");
      need(attempt.nominal.source.observer_exposure_second == f.bands[band].observer_exposure_second,
           "dark and photon exposure share optical exposure");
    }
  }
  near(*r.records[0].log_value, std::log(joint));
  near(*r.records[0].log_all_band_detection, std::log(selected));
  need(std::abs(joint - marginal0 * marginal1) > 1e-5L,
       "shared optical state creates cross-band dependence");
  f.records[0].measure = irred::detector::SelectionMeasure::selected_only;
  r = evaluate_optical_detector(f.input);
  need(r.status == S::ok, "selected-only joint event");
  near(*r.records[0].log_joint_record, std::log(joint));
  near(*r.records[0].log_value, std::log(joint / selected));
  f.observations[1] = {1, false, {}, {}};
  f.records[0].measure = irred::detector::SelectionMeasure::joint_detection_record;
  r = evaluate_optical_detector(f.input);
  joint = 0;
  for (unsigned s = 0; s < 2; ++s) {
    W n = f.full_photons() * (s ? .75L : .25L);
    joint += (s ? .75L : .25L) * poisson(n, 2) * std::exp(-n / 2);
  }
  need(r.status == S::ok && !r.attempts[0].nominal.rows[0].source.measured_adu,
       "nondetection retained");
  near(*r.records[0].log_value, std::log(joint));
  f.records[0].measure = irred::detector::SelectionMeasure::selected_only;
  r = evaluate_optical_detector(f.input);
  need(r.status == S::invalid_input && !r.records[0].log_value && r.conditional.size() == 2,
       "selected-only cannot drop a nondetection");
}
void normalization() {
  Fixture f;
  std::vector<std::array<irred::detector::Observation, 2>> data;
  std::vector<std::string> ids;
  std::vector<OpticalDetectorRecord> records;
  constexpr unsigned support = 24;
  data.reserve((support + 1) * (support + 1));
  ids.reserve(data.capacity()); records.reserve(data.capacity());
  for (unsigned a = 0; a <= support; ++a)
    for (unsigned b = 0; b <= support; ++b) {
      data.push_back({{{0, true, a, {}}, {0, true, b, {}}}});
      ids.push_back(std::to_string(a) + ":" + std::to_string(b));
    }
  for (std::size_t r = 0; r < data.size(); ++r) records.push_back({ids[r], data[r]});
  f.input.records = records;
  auto r = evaluate_optical_detector(f.input);
  need(r.status == S::ok && r.records.size() == data.size(), "coarse joint Poisson matrix");
  W native_sum = 0, exact_sum = 0;
  for (std::size_t i = 0; i < data.size(); ++i) {
    W ref = 0;
    for (unsigned s = 0; s < 2; ++s) {
      W n = f.full_photons() * (s ? .75L : .25L);
      ref += (s ? .75L : .25L) * poisson(n, *data[i][0].electron_count) *
             poisson(n / 2, *data[i][1].electron_count);
    }
    near(*r.records[i].log_value, std::log(ref));
    native_sum += std::exp(W(*r.records[i].log_value));
    exact_sum += ref;
  }
  need(1 - exact_sum < 1e-12L && exact_sum <= 1 + 1e-16L,
       "independent finite Poisson tail normalization");
  need(std::abs(native_sum - exact_sum) <= 2e-10L && native_sum > .999999999L,
       "joint mixture normalized including all count records");
}
void hostile() {
  Fixture f;
  auto base = evaluate_optical_detector(f.input);
  f.states[0].relative_mass = 2; f.states[1].relative_mass = 6;
  auto scaled = evaluate_optical_detector(f.input);
  near(*base.records[0].log_value, *scaled.records[0].log_value);
  std::swap(f.states[0], f.states[1]);
  auto permuted = evaluate_optical_detector(f.input);
  near(*base.records[0].log_value, *permuted.records[0].log_value);
  need(permuted.optical.normalized_state_mass[0] == .75L, "state permutation retained");
  std::swap(f.states[0], f.states[1]);
  f.detectors[1].quantum_efficiency = 2;
  auto invalid = evaluate_optical_detector(f.input);
  need(invalid.status == S::outside_domain && invalid.attempts.size() == 4 &&
       invalid.conditional.size() == 2 && !invalid.records[0].log_value,
       "required invalid detector state refuses whole mixture without dropping");
  need(invalid.attempts[0].nominal.status == S::ok && invalid.attempts[1].nominal.status == S::outside_domain,
       "successful other conditional attempts retained");
  f.detectors[1].quantum_efficiency = .5;
  f.detectors[0].background_expected_electrons =
      std::nextafter(64 - *base.optical.attempts[2].result.expected_photons.value, 0.0);
  f.observations[0] = {60, true, 64, {}};
  auto endpoint = evaluate_optical_detector(f.input);
  need(endpoint.status == S::outside_domain && !endpoint.records[0].log_value &&
       endpoint.attempts[2].nominal.status == S::ok &&
       endpoint.attempts[2].upper_photons->status == S::outside_domain,
       "required optical sensitivity endpoint outside detector domain is retained and refuses mixture");
  f.detectors[0].background_expected_electrons = 0;
  f.observations[0] = {1, true, 2, {}};
  f.high[1] = 1.1;
  invalid = evaluate_optical_detector(f.input);
  need(invalid.status == S::outside_domain && invalid.optical.attempts.size() == 4 &&
       invalid.records.empty(), "all optical state failures retained without detector mixture");
  f.high[1] = .75;
  f.input.conditional_law = ConditionalDetectorLaw::unspecified;
  need(evaluate_optical_detector(f.input).status == S::invalid_input, "conditional independence declaration required");
  f.input.conditional_law = ConditionalDetectorLaw::independent_poisson_and_read;
  f.observations[0].measure = irred::detector::SelectionMeasure::selected_only;
  invalid = evaluate_optical_detector(f.input);
  need(invalid.status == S::invalid_input && !invalid.records[0].log_value,
       "per-band selected denominators forbidden");
  f.observations[0].measure = irred::detector::SelectionMeasure::joint_detection_record;
  f.detectors[0].photon_law = irred::detector::PhotonLaw::unspecified;
  need(evaluate_optical_detector(f.input).status == S::outside_domain, "Poisson arrival declaration required");
  f.detectors[0].photon_law = irred::detector::PhotonLaw::poisson_arrivals;
  auto p = OpticalDetectorPolicy{}; p.maximum_payload_bytes = 1;
  need(evaluate_optical_detector(f.input, p).status == S::work_limit, "payload bound");
  p = {}; p.maximum_state_band_records = 3;
  need(evaluate_optical_detector(f.input, p).status == S::work_limit, "coarse cells bound");
  p = {}; p.detector_policy.maximum_poisson_terms = 1;
  invalid = evaluate_optical_detector(f.input, p);
  need(invalid.status == S::conditioning_budget_exceeded && invalid.attempts.size() == 4 &&
       !invalid.records[0].log_value, "parent term exhaustion and no-drop refusal");
  p = {}; p.absolute_log_allowance = 1e-30; p.scaled_log_allowance = 1e-30;
  need(evaluate_optical_detector(f.input, p).status == S::conditioning_budget_exceeded,
       "unmet final numerical budget refused");
  f.luminosity = {0, 0};
  f.observations = {{{1, false, {}, {}}, {1, false, {}, {}}}};
  auto zero = evaluate_optical_detector(f.input);
  need(zero.status == S::ok && zero.records[0].log_value && std::abs(*zero.records[0].log_value) < 2e-15 &&
       zero.records[0].zero_all_band_detection && !zero.attempts[0].lower_photons,
       "exact zero source nondetection law without fabricated positive event");
  f.observations = {{{0, true, 0, {}}, {0, true, 0, {}}}};
  f.records[0].measure = irred::detector::SelectionMeasure::selected_only;
  zero = evaluate_optical_detector(f.input);
  need(zero.status == S::ok && std::abs(*zero.records[0].log_value) < 2e-15 &&
       std::abs(*zero.records[0].log_all_band_detection) < 2e-15, "zero-noise exact threshold inclusion");
  f.observations[0].electron_count = 1;
  zero = evaluate_optical_detector(f.input);
  need(zero.status == S::ok && zero.records[0].zero_probability && !zero.records[0].log_value,
       "structural count-support zero retained");
  std::fesetround(FE_UPWARD);
  need(evaluate_optical_detector(f.input).status == S::conditioning_budget_exceeded,
       "arithmetic rounding contract");
  std::fesetround(FE_TONEAREST);
  f.observations = {{{20, false, {}, {}}, {20, false, {}, {}}}};
  f.detectors[0].read_noise_rms_electrons = f.detectors[1].read_noise_rms_electrons = 1;
  f.records[0].measure = irred::detector::SelectionMeasure::joint_detection_record;
  auto rare_event = evaluate_optical_detector(f.input);
  need(rare_event.status == S::conditioning_budget_exceeded &&
       !rare_event.records[0].log_all_band_detection && !rare_event.records[0].log_value &&
       rare_event.attempts[0].nominal.rows[0].status == S::ok,
       "published rare event diagnostic cannot bypass its gate on a joint nondetection record");
  f.detectors[0].read_noise_rms_electrons = f.detectors[1].read_noise_rms_electrons = 0;
  // An unresolved tiny spread fails the moment calculation, but not an explicitly
  // state-resolved likelihood using all actual conditional expectations.
  f.luminosity = {8e-24, 8e-24};
  f.low = {.5, .5}; f.high = {std::nextafter(.5, 1.0), std::nextafter(.5, 1.0)};
  f.records[0].measure = irred::detector::SelectionMeasure::joint_detection_record;
  f.observations = {{{1, true, 2, {}}, {1, true, 3, {}}}};
  need(evaluate_calibration(f.input.optical).status == S::conditioning_budget_exceeded,
       "unchanged tiny-spread moments refusal");
  need(evaluate_optical_detector(f.input).status == S::ok,
       "state likelihood requests no unrelated covariance approximation");
}
} // namespace
int main() {
  try {
    basic(); normalization(); hostile();
    std::cout << "PASS " << checks << " supplied optical joint detector controls\n";
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n'; return 1;
  }
}
