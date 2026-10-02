#include "irred/optical_detector.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>

namespace irred::photometry {
namespace {
using W = long double;
using S = numerics::Status;
constexpr W eps = std::numeric_limits<W>::epsilon();
constexpr W binary_eps = std::numeric_limits<double>::epsilon();
bool add(std::size_t &n, std::size_t x) {
  if (x > SIZE_MAX - n) return false;
  n += x;
  return true;
}
bool mul(std::size_t a, std::size_t b, std::size_t &n) {
  if (b && a > SIZE_MAX / b) return false;
  n = a * b;
  return true;
}
bool represent(W x, double &out) {
  out = static_cast<double>(x);
  return std::isfinite(x) && std::isfinite(out) && (x == 0 || std::isnormal(out));
}
struct LogLaw {
  S status = S::invalid_input;
  W record = 0, detection = 0, record_error = 0, detection_error = 0;
  bool zero_record = false, zero_detection = false;
};
LogLaw band_law(const detector::LikelihoodBatch &batch, std::size_t r,
                std::size_t records) {
  LogLaw out;
  if (batch.status != S::ok) { out.status = batch.status; return out; }
  if (batch.rows.size() != 2 * records) return out;
  const auto &data = batch.rows[r], &probe = batch.rows[records + r];
  if (data.status != S::ok) { out.status = data.status; return out; }
  if (probe.status != S::ok) { out.status = probe.status; return out; }
  if (!probe.detection_probability) return out;
  out.zero_record = data.zero_probability;
  if (!out.zero_record) {
    if (!data.log_value) return out;
    out.record = *data.log_value;
    out.record_error = data.numerical_error_estimate +
                       4 * binary_eps * (1 + std::abs(out.record));
  }
  const W detection = *probe.detection_probability;
  out.zero_detection = detection == 0;
  if (!out.zero_detection) {
    out.detection = std::log(detection);
    // The parent's nondetection log diagnostic is probability_error / below.
    // Its threshold-probability arithmetic is shared with the above tail.
    W absolute_error = batch.omitted_poisson_tail_estimate +
                       128 * eps * batch.poisson_terms;
    if (!probe.zero_probability) {
      if (!probe.log_value) return out;
      absolute_error = std::max(absolute_error,
          std::exp(W(*probe.log_value)) * probe.numerical_error_estimate);
    }
    out.detection_error = absolute_error / detection +
                          4 * binary_eps * (1 + std::abs(out.detection));
  }
  out.status = S::ok;
  return out;
}
LogLaw with_sensitivity(const OpticalDetectorAttempt &attempt, std::size_t r,
                        std::size_t records) {
  auto out = band_law(attempt.nominal, r, records);
  if (out.status != S::ok) return out;
  W record_sensitivity = 0, detection_sensitivity = 0;
  for (const auto *endpoint : {&attempt.lower_photons, &attempt.upper_photons}) {
    if (!*endpoint) continue;
    auto edge = band_law(**endpoint, r, records);
    if (edge.status != S::ok) { out.status = edge.status; return out; }
    if (edge.zero_record != out.zero_record || edge.zero_detection != out.zero_detection) {
      out.status = S::conditioning_budget_exceeded;
      return out;
    }
    if (!out.zero_record)
      record_sensitivity = std::max(record_sensitivity,
          std::abs(edge.record - out.record) + edge.record_error);
    if (!out.zero_detection)
      detection_sensitivity = std::max(detection_sensitivity,
          std::abs(edge.detection - out.detection) + edge.detection_error);
  }
  out.record_error += record_sensitivity;
  out.detection_error += detection_sensitivity;
  return out;
}
struct Mixture { W log = 0, error = 0; bool zero = true; };
Mixture mixture(std::span<const W> weights, const std::vector<LogLaw> &laws,
                std::size_t r, std::size_t records, bool event) {
  Mixture out;
  W maximum = -std::numeric_limits<W>::infinity();
  for (std::size_t s = 0; s < weights.size(); ++s) {
    const auto &law = laws[s * records + r];
    if (event ? law.zero_detection : law.zero_record) continue;
    const W value = std::log(weights[s]) + (event ? law.detection : law.record);
    maximum = std::max(maximum, value);
  }
  if (!std::isfinite(maximum)) return out;
  W total = 0, error = 0;
  for (std::size_t s = 0; s < weights.size(); ++s) {
    const auto &law = laws[s * records + r];
    if (event ? law.zero_detection : law.zero_record) continue;
    const W value = std::log(weights[s]) + (event ? law.detection : law.record);
    const W responsibility = std::exp(value - maximum);
    const W e = event ? law.detection_error : law.record_error;
    total += responsibility;
    error += responsibility * std::expm1(e);
  }
  out.zero = false;
  out.log = maximum + std::log(total);
  // Relative positive-mixture diagnostics become a log allocation; include
  // normalization, libm and wide accumulation ancestry explicitly.
  const W relative = error / total + (weights.size() + 16) * eps;
  out.error = relative < 1 ? -std::log1p(-relative)
                           : std::numeric_limits<W>::infinity();
  return out;
}
} // namespace

OpticalDetectorResult evaluate_optical_detector(const OpticalDetectorInput &in,
                                                OpticalDetectorPolicy p) {
  OpticalDetectorResult out;
  auto fail = [&](S status) {
    out.status = status;
    for (auto &r : out.records) {
      r.status = status;
      r.log_value.reset();
      r.log_joint_record.reset();
      r.log_all_band_detection.reset();
      r.zero_probability = r.zero_all_band_detection = false;
    }
    return std::move(out);
  };
  const auto ns = in.optical.states.size(), nb = in.optical.bands.size(),
             nr = in.records.size();
  if (std::numeric_limits<W>::digits < 64 || std::fegetround() != FE_TONEAREST)
    return fail(S::conditioning_budget_exceeded);
  if (!ns || !nb || !nr || in.detector_bands.size() != nb ||
      in.conditional_law != ConditionalDetectorLaw::independent_poisson_and_read ||
      in.detector_origin.empty() || in.observation_origin.empty() ||
      in.conditional_law_origin.empty() ||
      !std::isnormal(p.absolute_log_allowance) || p.absolute_log_allowance <= 0 ||
      !std::isnormal(p.scaled_log_allowance) || p.scaled_log_allowance <= 0)
    return fail(S::invalid_input);
  if (ns > 1024 || nb > 64 || nr > 1024 || ns > p.maximum_states ||
      nb > p.maximum_bands || nr > p.maximum_records ||
      p.maximum_payload_bytes > 1073741824 ||
      p.maximum_optical_payload_bytes > p.maximum_payload_bytes)
    return fail(S::work_limit);
  std::size_t attempts = 0, cells = 0, conditional = 0;
  if (!mul(ns, nb, attempts) || !mul(attempts, nr, cells) || !mul(ns, nr, conditional) ||
      cells > 65536 || cells > p.maximum_state_band_records)
    return fail(S::work_limit);
  for (std::size_t r = 0; r < nr; ++r) {
    if (in.records[r].id.empty() || in.records[r].bands.size() != nb) return out;
  }
  std::size_t bytes = sizeof(out);
  auto account = [&](std::size_t count, std::size_t size) {
    std::size_t amount = 0;
    return mul(count, size, amount) && add(bytes, amount);
  };
  const auto detector_bytes = detector::likelihood_payload_bound(2 * nr);
  if (!detector_bytes || !add(bytes, p.maximum_optical_payload_bytes) ||
      !account(attempts, sizeof(OpticalDetectorAttempt) + 3 * *detector_bytes) ||
      !account(conditional, sizeof(OpticalDetectorConditional) + sizeof(LogLaw)) ||
      !account(nr, sizeof(OpticalDetectorMixture) + 2 * sizeof(detector::Observation)))
    return fail(S::work_limit);
  for (auto origin : {in.detector_origin, in.observation_origin, in.conditional_law_origin})
    if (!add(bytes, origin.size()) || !add(bytes, 1)) return fail(S::work_limit);
  for (const auto &state : in.optical.states)
    if (!account(nb + nr, state.id.size()) || !account(nb + nr, 1)) return fail(S::work_limit);
  for (const auto &band : in.optical.bands)
    if (!account(ns, band.id.size()) || !account(ns, 1)) return fail(S::work_limit);
  for (const auto &record : in.records)
    if (!account(ns + 1, record.id.size()) || !account(ns + 1, 1)) return fail(S::work_limit);
  if (bytes > p.maximum_payload_bytes) return fail(S::work_limit);
  for (std::size_t r = 0; r < nr; ++r)
    for (std::size_t earlier = 0; earlier < r; ++earlier)
      if (in.records[r].id == in.records[earlier].id) return out;
  try {
    CalibrationPolicy optical_policy;
    optical_policy.requested_outputs = transmitted_photons;
    optical_policy.aggregation = CalibrationAggregation::state_resolved_only;
    optical_policy.maximum_states = p.maximum_states;
    optical_policy.maximum_bands = p.maximum_bands;
    optical_policy.maximum_total_knots = p.maximum_total_optical_knots;
    optical_policy.maximum_output_bytes = p.maximum_optical_payload_bytes;
    out.optical = evaluate_calibration(in.optical, optical_policy);
    if (out.optical.status != S::ok) return fail(out.optical.status);
    out.detector_origin = in.detector_origin;
    out.observation_origin = in.observation_origin;
    out.conditional_law_origin = in.conditional_law_origin;
    out.conditional_law = in.conditional_law;
    out.attempts.reserve(attempts);
    out.conditional.reserve(conditional);
    out.records.reserve(nr);
    std::vector<detector::Observation> observations(2 * nr);
    for (std::size_t s = 0; s < ns; ++s)
      for (std::size_t b = 0; b < nb; ++b) {
        for (std::size_t r = 0; r < nr; ++r) {
          observations[r] = in.records[r].bands[b];
          observations[nr + r] = {observations[r].detection_threshold_adu, false, {}, {},
                                  detector::SelectionMeasure::joint_detection_record};
        }
        const auto &band = in.detector_bands[b];
        const auto photons = *out.optical.attempts[s * nb + b].result.expected_photons.value;
        detector::Input source{band.photon_law, photons, band.quantum_efficiency,
                              band.background_expected_electrons, band.dark_electrons_per_second,
                              in.optical.bands[b].observer_exposure_second,
                              band.read_noise_rms_electrons, band.gain_electrons_per_adu, band.bias_adu};
        OpticalDetectorAttempt attempt;
        attempt.state_id = in.optical.states[s].id;
        attempt.band_id = in.optical.bands[b].id;
        attempt.nominal = detector::likelihood(source, observations, p.detector_policy);
        out.poisson_terms += attempt.nominal.poisson_terms;
        if (photons > 0) {
          const W allocation = calibration_sampled_relative_sensitivity * photons;
          double lower = std::nextafter(double(W(photons) - allocation), 0.0);
          double upper = std::nextafter(double(W(photons) + allocation),
                                       std::numeric_limits<double>::infinity());
          for (unsigned edge = 0; edge < 2; ++edge) {
            const double n = edge ? upper : lower;
            detector::LikelihoodBatch endpoint;
            source.expected_transmitted_photons = n;
            if (!std::isnormal(n) || n <= 0) {
              endpoint.source = source;
              endpoint.status = S::outside_domain;
            } else endpoint = detector::likelihood(source, observations, p.detector_policy);
            out.poisson_terms += endpoint.poisson_terms;
            (edge ? attempt.upper_photons : attempt.lower_photons) = std::move(endpoint);
          }
        }
        out.attempts.push_back(std::move(attempt));
      }
    std::vector<LogLaw> laws(conditional);
    S aggregate_status = S::ok;
    for (std::size_t s = 0; s < ns; ++s)
      for (std::size_t r = 0; r < nr; ++r) {
        auto &law = laws[s * nr + r];
        law.status = S::ok;
        for (std::size_t b = 0; b < nb; ++b) {
          auto value = with_sensitivity(out.attempts[s * nb + b], r, nr);
          if (value.status != S::ok) { law.status = value.status; continue; }
          law.record += value.record;
          law.detection += value.detection;
          law.record_error += value.record_error;
          law.detection_error += value.detection_error;
          law.zero_record = law.zero_record || value.zero_record;
          law.zero_detection = law.zero_detection || value.zero_detection;
          const auto &record = in.records[r];
          if (record.bands[b].measure != detector::SelectionMeasure::joint_detection_record ||
              (record.measure == detector::SelectionMeasure::selected_only && !record.bands[b].detected))
            law.status = S::invalid_input;
        }
        if (in.records[r].measure != detector::SelectionMeasure::joint_detection_record &&
            in.records[r].measure != detector::SelectionMeasure::selected_only)
          law.status = S::invalid_input;
        OpticalDetectorConditional c;
        c.state_id = in.optical.states[s].id;
        c.record_id = in.records[r].id;
        c.status = law.status;
        if (c.status == S::ok) {
          double log = 0;
          if (!law.zero_record) {
            if (!represent(law.record, log)) c.status = S::outside_domain;
            else { c.log_record = log; law.record_error += std::abs(W(log) - law.record); }
          }
          if (!law.zero_detection) {
            if (!represent(law.detection, log)) c.status = S::outside_domain;
            else { c.log_all_band_detection = log; law.detection_error += std::abs(W(log) - law.detection); }
          }
          if (!represent(law.record_error, c.record_log_error_estimate) ||
              !represent(law.detection_error, c.detection_log_error_estimate)) c.status = S::outside_domain;
          c.zero_record = law.zero_record;
          c.zero_all_band_detection = law.zero_detection;
        }
        law.status = c.status;
        out.conditional.push_back(std::move(c));
      }
    for (std::size_t r = 0; r < nr; ++r) {
      OpticalDetectorMixture record;
      record.record_id = in.records[r].id;
      record.measure = in.records[r].measure;
      record.status = S::ok;
      for (std::size_t s = 0; s < ns; ++s)
        if (laws[s * nr + r].status != S::ok) record.status = laws[s * nr + r].status;
      if (record.status == S::ok) {
        const auto joint = mixture(out.optical.normalized_state_mass, laws, r, nr, false),
                   event = mixture(out.optical.normalized_state_mass, laws, r, nr, true);
        const bool selected = record.measure == detector::SelectionMeasure::selected_only;
        double log = 0;
        if (!joint.zero) {
          if (!represent(joint.log, log)) record.status = S::outside_domain;
          else {
            const W error = joint.error + std::abs(W(log) - joint.log);
            if (!std::isfinite(error) ||
                error > p.absolute_log_allowance + p.scaled_log_allowance * (1 + std::abs(joint.log)))
              record.status = S::conditioning_budget_exceeded;
            else if (!represent(error, record.joint_record_log_error_estimate))
              record.status = S::outside_domain;
            else record.log_joint_record = log;
          }
        }
        if (!event.zero) {
          if (!represent(event.log, log)) record.status = S::outside_domain;
          else {
            const W error = event.error + std::abs(W(log) - event.log);
            if (!std::isfinite(error) ||
                error > p.absolute_log_allowance + p.scaled_log_allowance * (1 + std::abs(event.log)))
              record.status = S::conditioning_budget_exceeded;
            else if (!represent(error, record.detection_log_error_estimate))
              record.status = S::outside_domain;
            else record.log_all_band_detection = log;
          }
        }
        record.zero_all_band_detection = event.zero;
        if (selected && event.zero) record.status = S::outside_domain;
        if (!joint.zero && record.status == S::ok) {
          const W value = joint.log - (selected ? event.log : 0);
          W error = joint.error + (selected ? event.error : 0);
          if (!represent(value, log)) record.status = S::outside_domain;
          else {
            error += std::abs(W(log) - value);
            if (!std::isfinite(error) ||
                error > p.absolute_log_allowance + p.scaled_log_allowance * (1 + std::abs(value)))
              record.status = S::conditioning_budget_exceeded;
            else if (!represent(error, record.numerical_error_estimate)) record.status = S::outside_domain;
            else record.log_value = log;
          }
        }
        if (joint.zero && record.status == S::ok) record.zero_probability = true;
      }
      if (record.status != S::ok) aggregate_status = record.status;
      out.records.push_back(std::move(record));
    }
    out.status = aggregate_status;
    return out;
  } catch (const std::bad_alloc &) { return fail(S::work_limit); }
    catch (const std::length_error &) { return fail(S::work_limit); }
}
} // namespace irred::photometry
