#include "irred/detector_selection.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred::detector;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void need(bool v, const char *why) {
  ++checks;
  if (!v)
    throw std::runtime_error(why);
}
void near(long double a, long double b, long double budget = 2e-12L) {
  need(std::abs(a - b) <= budget * (1 + std::abs(b)),
       "frozen analytic allocation");
}
Input input() { return {PhotonLaw::poisson_arrivals, 4, .5, 0, 0, 1, 0, 1, 0}; }
SelectionPolicy policy() {
  return {65536, 64 * 1024 * 1024, 256, 5e-13, 2e-11};
}
void known_answers() {
  // Three factual known-answer rows,
  // Random1239545ff6413f258be2f04c1d319d99aaef7521150/tests/kat_vectors.
  need(random_words(0, {0, 0}) ==
           std::array<std::uint32_t, 4>{0x6627e8d5, 0xe169c58d, 0xbc57ac4c,
                                        0x9b00dbd8},
       "Philox zero KAT");
  need(random_words(UINT64_MAX, {UINT64_MAX, UINT64_MAX}) ==
           std::array<std::uint32_t, 4>{0x408f276d, 0x41c83b0e, 0xa20bc7c6,
                                        0x6d5451fd},
       "Philox max KAT");
  need(random_words(0x299f31d0a4093822ULL,
                    {0x0370734413198a2eULL, 0x85a308d3243f6a88ULL}) ==
           std::array<std::uint32_t, 4>{0xd16cfe09, 0x94fdcceb, 0x5001e420,
                                        0x24126ea1},
       "Philox mixed KAT");
}
void exact() {
  auto x = input();
  auto m = moments(x);
  need(m.status == S::ok, "moments admitted");
  near(*m.signal_electrons, 2);
  near(*m.poisson_mean_electrons, 2);
  near(*m.mean_adu, 2);
  near(*m.variance_adu, 2);
  x.background_expected_electrons = 1;
  x.dark_electrons_per_second = .25;
  x.observer_exposure_second = 4;
  x.read_noise_rms_electrons = 2;
  x.gain_electrons_per_adu = 2;
  x.bias_adu = 3;
  m = moments(x);
  near(*m.poisson_mean_electrons, 4);
  near(*m.mean_adu, 5);
  near(*m.variance_adu, 2);
  x = input();
  const Observation observations[] = {
      {3, false, {}, {}, SelectionMeasure::joint_detection_record},
      {3, true, 3, {}, SelectionMeasure::joint_detection_record},
      {3, true, 3, {}, SelectionMeasure::selected_only},
      {3, true, 2, {}, SelectionMeasure::joint_detection_record}};
  auto result = likelihood(x, observations, policy());
  need(result.status == S::ok && result.rows.size() == 4, "selection batch");
  const long double nondetected = 5 * std::exp(-2.L),
                    mass = 4.L / 3 * std::exp(-2.L), detected = 1 - nondetected;
  for (unsigned j = 0; j < 3; ++j)
    need(result.rows[j].status == S::ok, "exact Poisson likelihood admitted");
  near(*result.rows[0].log_value, std::log(nondetected));
  near(*result.rows[0].detection_probability, detected);
  near(*result.rows[1].log_value, std::log(mass));
  near(*result.rows[2].log_value, std::log(mass / detected));
  need(result.rows[3].zero_probability && !result.rows[3].log_value,
       "inconsistent detected record is outside support");
  x.expected_transmitted_photons = 0;
  x.read_noise_rms_electrons = 1;
  x.gain_electrons_per_adu = 2;
  x.bias_adu = 3;
  const Observation gaussian[] = {
      {3, false, {}, {}, SelectionMeasure::joint_detection_record},
      {3, true, {}, 3.5, SelectionMeasure::joint_detection_record},
      {3, true, {}, 3.5, SelectionMeasure::selected_only}};
  result = likelihood(x, gaussian, policy());
  for (const auto &r : result.rows)
    need(r.status == S::ok, "zero-Poisson Gaussian admitted");
  near(*result.rows[0].log_value, -std::log(2.L));
  near(*result.rows[1].log_value,
       std::log(2.L) - .5L -
           .5L * std::log(2 * std::numbers::pi_v<long double>));
  near(*result.rows[2].log_value, *result.rows[1].log_value + std::log(2.L));
  x.read_noise_rms_electrons = 0;
  const Observation deterministic[] = {
      {4, false, {}, {}, SelectionMeasure::joint_detection_record},
      {3, false, {}, {}, SelectionMeasure::joint_detection_record},
      {4, true, 0, {}, SelectionMeasure::joint_detection_record}};
  result = likelihood(x, deterministic, policy());
  need(result.rows[0].status == S::ok && *result.rows[0].log_value == 0,
       "deterministic censor atom");
  need(result.rows[1].zero_probability && result.rows[2].zero_probability,
       "structural atom zeros explicit");
}
void sampling() {
  constexpr unsigned count = 65536;
  auto x = input();
  x.expected_transmitted_photons = 12;
  x.quantum_efficiency = .75;
  x.background_expected_electrons = 1;
  x.read_noise_rms_electrons = 2;
  x.gain_electrons_per_adu = 2;
  x.bias_adu = 3;
  std::vector<Request> requests;
  requests.reserve(count);
  for (unsigned i = 0; i < count; ++i)
    requests.push_back({x, {19, i}});
  const auto bytes = output_payload_bound(count);
  need(bytes.has_value(), "draw payload bound");
  auto batch = simulate(requests, 0x12345678, {count, *bytes});
  need(batch.status == S::ok && batch.rows.size() == count,
       "bounded ensemble admitted");
  long double sum = 0, sum2 = 0, psum = 0, rsum = 0, rsum2 = 0, cross = 0;
  unsigned below = 0;
  for (const auto &r : batch.rows) {
    need(r.status == S::ok && r.measured_adu && r.poisson_electrons &&
             r.read_electrons,
         "all addresses retained");
    long double y = *r.measured_adu, p = *r.poisson_electrons,
                read = *r.read_electrons;
    sum += y;
    sum2 += y * y;
    psum += p;
    rsum += read;
    rsum2 += read * read;
    cross += (p - 10) * read;
    below += p < 10;
  }
  const long double mean = sum / count, variance = sum2 / count - mean * mean;
  // Predeclared6-sigma finite-ensemble controls. E(Y)=8, Var(Y)=3.5.
  near(mean, 8,
       6 * std::sqrt(3.5L / count) / 9); // near scales by1+abs(reference)
  const long double fourth = 3 * 3.5L * 3.5L + 10.L / 16;
  need(std::abs(variance - 3.5L) <=
           6 * std::sqrt((fourth - 3.5L * 3.5L) / count),
       "variance analytic fourth-moment control");
  need(std::abs(psum / count - 10) <= 6 * std::sqrt(10.L / count),
       "Poisson mean control");
  need(std::abs(rsum / count) <= 6 * 2 / std::sqrt((long double)count),
       "Gaussian mean control");
  need(std::abs(rsum2 / count - 4) <= 6 * std::sqrt(32.L / count),
       "Gaussian second-moment control");
  need(std::abs(cross / count) <= 6 * std::sqrt(40.L / count),
       "stream-component empirical cross-moment");
  long double term = std::exp(-10.L), cdf = term;
  for (unsigned n = 1; n < 10; ++n) {
    term *= 10.L / n;
    cdf += term;
  }
  need(std::abs((long double)below / count - cdf) <=
           6 * std::sqrt(cdf * (1 - cdf) / count) + std::ldexp(1.L, -32),
       "Poisson CDF and finite uniform control");
  auto replay =
      simulate(std::span(requests.data(), 3), 0x12345678, {3, 1024 * 1024});
  for (unsigned j = 0; j < 3; ++j)
    need(replay.rows[j].words == batch.rows[j].words &&
             replay.rows[j].measured_adu == batch.rows[j].measured_adu,
         "address replay independent of batch size");
  std::swap(requests[0], requests[2]);
  replay =
      simulate(std::span(requests.data(), 3), 0x12345678, {3, 1024 * 1024});
  need(replay.rows[0].words == batch.rows[2].words, "sharding/order replay");
  requests[1].address = requests[0].address;
  need(
      simulate(std::span(requests.data(), 3), 0, {3, 1024 * 1024}).rows.empty(),
      "duplicate addresses rejected");
}
void invalid() {
  auto x = input();
  x.photon_law = PhotonLaw::unspecified;
  need(moments(x).status != S::ok, "photon generating law explicit");
  x = input();
  x.expected_transmitted_photons = 129;
  need(moments(x).status == S::outside_domain, "large Poisson unsupported");
  x = input();
  x.expected_transmitted_photons = std::numeric_limits<double>::denorm_min();
  need(!moments(x).mean_adu, "subnormal source refused");
  x = input();
  x.gain_electrons_per_adu = 1e308;
  need(moments(x).status == S::outside_domain,
       "positive variance cannot fabricate zero");
  x = input();
  x.expected_transmitted_photons = 2;
  x.gain_electrons_per_adu = 1 + std::ldexp(1., -52);
  x.bias_adu = -(1 - std::ldexp(1., -52));
  need(moments(x).status != S::ok,
       "unresolved positive bias cancellation cannot fabricate zero");
  x = input();
  x.read_noise_rms_electrons = 1;
  const Observation bad[] = {
      {0, true, 1, {}},
      {0, false, 1, {}},
      {0, false, {}, {}, SelectionMeasure::selected_only}};
  auto r = likelihood(x, bad, policy());
  for (const auto &row : r.rows)
    need(row.status == S::invalid_input, "measure/payload mismatch refused");
  auto p = policy();
  p.maximum_poisson_terms = 1;
  const Observation o{
      3, false, {}, {}, SelectionMeasure::joint_detection_record};
  r = likelihood(input(), std::span(&o, 1), p);
  need(!r.rows[0].log_value, "unresolved omitted tail fails");
  p = policy();
  p.maximum_payload_bytes = *likelihood_payload_bound(1) - 1;
  need(likelihood(input(), std::span(&o, 1), p).status == S::work_limit,
       "likelihood payload admission");
  const Request q{input(), {0, 0}};
  need(
      simulate(std::span(&q, 1), 0, {1, *output_payload_bound(1) - 1}).status ==
          S::work_limit,
      "simulation payload admission");
  need(!output_payload_bound(SIZE_MAX) && !likelihood_payload_bound(SIZE_MAX),
       "payload overflow");
  const int rounding = std::fegetround();
  std::fesetround(FE_UPWARD);
  need(moments(input()).status != S::ok, "arithmetic contract");
  std::fesetround(rounding);
}
void censored_recovery() {
  constexpr unsigned n = 8192;
  const auto x = input();
  std::vector<Request> requests;
  for (unsigned j = 0; j < n; ++j)
    requests.push_back({x, {55, j}});
  const auto draws = simulate(requests, 42, {n, *output_payload_bound(n)});
  std::vector<Observation> observations;
  unsigned detections = 0;
  long double detected_sum = 0;
  for (const auto &row : draws.rows) {
    need(row.status == S::ok, "synthetic censor generating law");
    const auto count = *row.poisson_electrons;
    observations.push_back(
        {3,
         count >= 3,
         count >= 3 ? std::optional<std::uint32_t>(count) : std::nullopt,
         {},
         SelectionMeasure::joint_detection_record});
    if (count >= 3) {
      ++detections;
      detected_sum += count;
    }
  }
  need(detections > 0 && detections < n, "non-detections retained in recovery");
  // External consumer solves this toy censored Poisson score independently.
  // F_2(lambda)=exp(-lambda)(1+lambda+lambda^2/2).
  auto score = [&](long double l) {
    return detected_sum / l - detections -
           (n - detections) * l * l / (2 + 2 * l + l * l);
  };
  long double low = .1L, high = 10;
  for (unsigned j = 0; j < 80; ++j) {
    const auto mid = (low + high) / 2;
    if (score(mid) > 0)
      low = mid;
    else
      high = mid;
  }
  const long double recovered = (low + high) / 2;
  // Frozen conservative six-sigma allocation for this finite seeded sample.
  need(std::abs(recovered - 2) < .15,
       "censored joint recovery of synthetic mean");
  need(detected_sum / detections > 3,
       "dropping non-detections biases raw mean");
  auto objective = [&](long double l) {
    auto model = x;
    model.expected_transmitted_photons = 2 * double(l);
    const auto evaluated = likelihood(model, observations, policy());
    need(evaluated.rows.size() == n, "all censored records evaluated");
    long double total = 0;
    for (const auto &row : evaluated.rows) {
      need(row.status == S::ok && row.log_value,
           "joint recovery likelihood finite");
      total += *row.log_value;
    }
    return total;
  };
  const auto best = objective(recovered);
  need(best > objective(recovered - .01L) && best > objective(recovered + .01L),
       "native full-record objective agrees with external score maximum");
}
} // namespace
int main() {
  try {
    known_answers();
    exact();
    sampling();
    invalid();
    censored_recovery();
    std::cout << "PASS " << checks << " detector/selection controls\n";
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
