#include "irred/detector_selection.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
namespace irred::detector {
namespace {
using W = long double;
using S = numerics::Status;
constexpr W eps = std::numeric_limits<W>::epsilon();
bool normal(double x) {
  return std::isfinite(x) && (x == 0 || std::fpclassify(x) == FP_NORMAL);
}
bool accepted(W x, W error, const SelectionPolicy &p) {
  return x >= 0 && std::isfinite(x) &&
         error <= W(p.absolute_probability_allowance) +
                      p.relative_probability_allowance * x;
}
bool represent(W x, double &d) {
  d = static_cast<double>(x);
  return std::isfinite(d) &&
         (x == 0 || (d != 0 && std::fpclassify(d) == FP_NORMAL));
}
struct Mapped { W value, error; };
// TwoSum subtraction and FMA product residual retain a lost threshold edge.
// The returned residual/arithmetic diagnostic is not a rigorous libm bound.
Mapped electron_coordinate(W adu, W bias, W gain) {
  const W difference = adu - bias, bp = difference - adu;
  const W low = (adu - (difference - bp)) + (-bias - bp);
  const W high = difference * gain;
  const W product_low = std::fma(difference, gain, -high);
  const W low_product = low * gain;
  const W correction = product_low + low_product, value = high + correction;
  const W error = std::abs(correction - (value - high)) +
                  8 * eps * (std::abs(product_low) + std::abs(low_product));
  return {value, error};
}
} // namespace
std::optional<std::size_t> likelihood_payload_bound(std::size_t n) noexcept {
  if (n > (SIZE_MAX - sizeof(LikelihoodBatch)) / sizeof(LikelihoodRow))
    return {};
  return sizeof(LikelihoodBatch) + n * sizeof(LikelihoodRow);
}
LikelihoodBatch likelihood(const Input &source,
                           std::span<const Observation> observations,
                           SelectionPolicy p) {
  LikelihoodBatch out;
  out.source = source;
  const auto m = moments(source);
  if (m.status != S::ok) {
    out.status = m.status;
    return out;
  }
  if (p.maximum_rows > 65536 || observations.size() > 65536 ||
      !p.maximum_poisson_terms || p.maximum_poisson_terms > 256 ||
      !normal(p.absolute_probability_allowance) ||
      !normal(p.relative_probability_allowance) ||
      p.absolute_probability_allowance <= 0 ||
      p.relative_probability_allowance <= 0)
    return out;
  const auto b = likelihood_payload_bound(observations.size());
  if (!b || observations.size() > p.maximum_rows ||
      *b > p.maximum_payload_bytes) {
    out.status = S::work_limit;
    return out;
  }
  const W lambda = *m.poisson_mean_electrons,
          sd = source.read_noise_rms_electrons,
          g = source.gain_electrons_per_adu, bias = source.bias_adu;
  std::array<W, 256> weights{};
  weights[0] = std::exp(-lambda);
  W tail = 1;
  for (std::size_t n = 0; n < p.maximum_poisson_terms; ++n) {
    if (n)
      weights[n] = weights[n - 1] * lambda / n;
    out.poisson_terms = n + 1;
    if (lambda == 0) {
      tail = 0;
      break;
    }
    if (W(n + 2) > lambda) {
      tail = weights[n] * lambda / (n + 1) / (1 - lambda / (n + 2));
      if (tail < 1e-28L)
        break;
    }
  }
  if (!represent(tail, out.omitted_poisson_tail_estimate)) {
    out.status = S::outside_domain;
    return out;
  }
  out.status = S::ok;
  out.rows.reserve(observations.size());
  for (const auto &obs : observations) {
    out.rows.emplace_back();
    auto &r = out.rows.back();
    r.source = obs;
    if (!normal(obs.detection_threshold_adu) ||
        (obs.measure != SelectionMeasure::joint_detection_record &&
         obs.measure != SelectionMeasure::selected_only) ||
        (!obs.detected && (obs.electron_count || obs.measured_adu ||
                           obs.measure == SelectionMeasure::selected_only)) ||
        (obs.detected &&
         ((sd == 0 && (!obs.electron_count || obs.measured_adu)) ||
          (sd > 0 && (!obs.measured_adu || obs.electron_count)))) ||
        (obs.measured_adu && !normal(*obs.measured_adu)))
      continue;
    const auto mapped_threshold = electron_coordinate(obs.detection_threshold_adu,bias,g);
    const W threshold = mapped_threshold.value;
    const auto mapped_observed = obs.measured_adu
      ? electron_coordinate(*obs.measured_adu,bias,g) : Mapped{0,0};
    if (!std::isfinite(threshold) || !std::isfinite(mapped_observed.value)) {
      r.status=S::overflow; continue;
    }
    W below = 0, above = 0, density = 0;
    W mapping_density_error=0;
    bool unresolved=false;
    // Sum both tails directly: 1-CDF would lose rare detections.
    for (std::size_t n = 0; n < out.poisson_terms; ++n) {
      if (sd == 0) {
        if (mapped_threshold.error>0 &&
            std::abs(W(n)-threshold)<=mapped_threshold.error) {
          unresolved=true; break;
        }
        if (W(n) < threshold)
          below += weights[n];
        else
          above += weights[n];
      } else {
        const W t = (threshold - W(n)) / (sd * std::sqrt(2.L));
        below += weights[n] * std::erfc(-t) / 2;
        above += weights[n] * std::erfc(t) / 2;
        if (obs.detected) {
          const W z = (mapped_observed.value-W(n)) / sd;
          const W term=weights[n] * std::exp(-z*z/2)*g /
                     (sd*std::sqrt(2*std::numbers::pi_v<W>));
          const W dz=mapped_observed.error/sd;
          const W log_error=std::abs(z)*dz+dz*dz/2;
          if (!std::isfinite(log_error) || log_error>.1L) { unresolved=true; break; }
          density+=term;
          mapping_density_error+=term*std::expm1(log_error);
        }
      }
    }
    if (sd==0 && obs.electron_count && mapped_threshold.error>0 &&
        std::abs(W(*obs.electron_count)-threshold)<=mapped_threshold.error)
      unresolved=true;
    if (unresolved) { r.status=S::conditioning_budget_exceeded; continue; }
    const W mapping_probability_error=sd>0
      ? mapped_threshold.error/(sd*std::sqrt(2*std::numbers::pi_v<W>)) : 0;
    const W probability_error = tail + 128 * eps * out.poisson_terms + mapping_probability_error;
    W chosen = obs.detected ? density : below;
    W chosen_error = probability_error;
    bool structural_zero = false;
    if (sd == 0) {
      if (!obs.detected && threshold <= 0)
        structural_zero = true;
      if (obs.detected) {
        const W n = *obs.electron_count;
        if (n < threshold)
          structural_zero = true;
        if (lambda == 0 && n > 0)
          structural_zero = true;
        if (!structural_zero) {
          const W logp =
              lambda == 0 ? 0
                          : -lambda + n * std::log(lambda) - std::lgamma(n + 1);
          chosen = std::exp(logp);
          chosen_error = std::abs(chosen) * 128 * eps * (1 + std::abs(logp));
        }
      }
    } else {
      if (obs.detected && *obs.measured_adu < obs.detection_threshold_adu)
        structural_zero = true;
      if (obs.detected)
        chosen_error = tail * g / (sd * std::sqrt(2 * std::numbers::pi_v<W>)) +
                       128 * eps * out.poisson_terms * std::abs(density) + mapping_density_error;
    }
    if (!accepted(above, probability_error, p) ||
        !accepted(chosen, chosen_error, p)) {
      r.status = S::conditioning_budget_exceeded;
      continue;
    }
    double detection = 0;
    if (above == 0 && !(sd == 0 && lambda == 0 && threshold > 0)) {
      r.status = S::outside_domain;
      continue;
    }
    if (!represent(above, detection)) {
      r.status = S::outside_domain;
      continue;
    }
    r.detection_probability = detection;
    if (obs.measure == SelectionMeasure::selected_only && above == 0) {
      r.status = S::outside_domain;
      continue;
    }
    if (structural_zero) {
      r.status = S::ok;
      r.zero_probability = true;
      continue;
    }
    if (chosen <= 0) {
      r.status = S::outside_domain;
      continue;
    }
    W logp = std::log(chosen), error = chosen_error / chosen;
    if (obs.measure == SelectionMeasure::selected_only) {
      logp -= std::log(above);
      error += probability_error / above;
    }
    if (!std::isfinite(logp) ||
        error > p.relative_probability_allowance +
                    p.absolute_probability_allowance / (1 + std::abs(logp))) {
      r.status = S::conditioning_budget_exceeded;
      continue;
    }
    double result = 0, diagnostic = 0;
    if (!represent(logp, result) || !represent(error, diagnostic)) {
      r.status = S::outside_domain;
      continue;
    }
    r.status = S::ok;
    r.log_value = result;
    r.numerical_error_estimate = diagnostic;
  }
  return out;
}
} // namespace irred::detector
