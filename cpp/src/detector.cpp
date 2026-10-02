#include "irred/detector.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::detector {
namespace {
using W = long double;
using S = numerics::Status;
bool arithmetic() {
  return std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384 &&
         std::fegetround() == FE_TONEAREST;
}
S store(W x, std::optional<double> &v) {
  const double d = static_cast<double>(x);
  if (!std::isfinite(x) || !std::isfinite(d))
    return S::overflow;
  if (x != 0 && (d == 0 || std::fpclassify(d) != FP_NORMAL))
    return S::outside_domain;
  v = d;
  return S::ok;
}
S biased_store(W numerator, W low, W gain, W bias, std::optional<double> &out) {
  const W q = numerator / gain, value = q + bias;
  const W product = -bias * gain;
  const bool exact_zero =
      low == 0 && numerator == product && std::fma(-bias, gain, -product) == 0;
  const W estimate = std::abs(low / gain) +
                     8 * std::numeric_limits<W>::epsilon() *
                         (std::abs(q) + std::abs(bias) + std::abs(value));
  if (value == 0) {
    if (!exact_zero)
      return S::conditioning_budget_exceeded;
  } else if (estimate > 2e-12L * std::abs(value))
    return S::conditioning_budget_exceeded;
  return store(value, out);
}
S admit(const Input &p) {
  for (double x : {p.expected_transmitted_photons, p.quantum_efficiency,
                   p.background_expected_electrons, p.dark_electrons_per_second,
                   p.observer_exposure_second, p.read_noise_rms_electrons,
                   p.gain_electrons_per_adu, p.bias_adu}) {
    if (!std::isfinite(x))
      return S::nonfinite_input;
    if (x != 0 && std::fpclassify(x) != FP_NORMAL)
      return S::outside_domain;
  }
  if (p.photon_law != PhotonLaw::poisson_arrivals ||
      p.expected_transmitted_photons < 0 || p.quantum_efficiency < 0 ||
      p.quantum_efficiency > 1 || p.background_expected_electrons < 0 ||
      p.dark_electrons_per_second < 0 || p.observer_exposure_second < 0 ||
      p.read_noise_rms_electrons < 0 || p.read_noise_rms_electrons > 64 ||
      p.gain_electrons_per_adu <= 0)
    return S::outside_domain;
  return S::ok;
}
} // namespace
Moments moments(const Input &p) noexcept {
  Moments out;
  if (!arithmetic())
    return out;
  out.status = admit(p);
  if (out.status != S::ok)
    return out;
  const W signal = W(p.quantum_efficiency) * p.expected_transmitted_photons;
  const W lambda = signal + p.background_expected_electrons +
                   W(p.dark_electrons_per_second) * p.observer_exposure_second;
  if (!std::isfinite(lambda) || lambda > 64) {
    out.status = S::outside_domain;
    return out;
  }
  const W g = p.gain_electrons_per_adu, sd = p.read_noise_rms_electrons;
  out.status = store(signal, out.signal_electrons);
  if (out.status == S::ok)
    out.status = store(lambda, out.poisson_mean_electrons);
  if (out.status == S::ok)
    out.status = biased_store(*out.poisson_mean_electrons, 0, g, p.bias_adu,
                              out.mean_adu);
  if (out.status == S::ok)
    out.status = store((W(*out.poisson_mean_electrons) + sd * sd) / (g * g),
                       out.variance_adu);
  if (out.status != S::ok) {
    out.signal_electrons.reset();
    out.poisson_mean_electrons.reset();
    out.mean_adu.reset();
    out.variance_adu.reset();
  }
  return out;
}
std::array<std::uint32_t, 4> random_words(std::uint64_t seed,
                                          Address a) noexcept {
  return random::words(seed, a);
}
std::optional<std::size_t> output_payload_bound(std::size_t n) noexcept {
  if (n > (SIZE_MAX - sizeof(Batch)) / (sizeof(Draw) + sizeof(Address)))
    return {};
  return sizeof(Batch) + n * (sizeof(Draw) + sizeof(Address));
}
Batch simulate(std::span<const Request> inputs, std::uint64_t seed, Policy p) {
  Batch out;
  out.seed = seed;
  if (!arithmetic() || p.maximum_rows > 65536 || inputs.size() > 65536)
    return out;
  const auto b = output_payload_bound(inputs.size());
  if (!b || inputs.size() > p.maximum_rows || *b > p.maximum_payload_bytes) {
    out.status = S::work_limit;
    return out;
  }
  std::vector<Address> addresses;
  addresses.reserve(inputs.size());
  for (const auto &x : inputs)
    addresses.push_back(x.address);
  auto less = [](Address a, Address b) {
    return a.stream < b.stream || (a.stream == b.stream && a.sample < b.sample);
  };
  std::sort(addresses.begin(), addresses.end(), less);
  for (std::size_t j = 1; j < addresses.size(); ++j)
    if (!less(addresses[j - 1], addresses[j]))
      return out;
  out.rows.reserve(inputs.size());
  out.status = S::ok;
  for (const auto &x : inputs) {
    out.rows.emplace_back();
    auto &d = out.rows.back();
    d.source = x.source;
    d.address = x.address;
    d.words = random_words(seed, x.address);
    const auto m = moments(x.source);
    d.status = m.status;
    if (m.status != S::ok)
      continue;
    const W lambda = *m.poisson_mean_electrons;
    W weight = std::exp(-lambda), cdf = weight;
    const W u = random::halfbin(d.words[0]);
    unsigned n = 0;
    while (u > cdf && n < 255) {
      ++n;
      weight *= lambda / n;
      cdf += weight;
    }
    if (u > cdf) {
      d.status = S::work_limit;
      continue;
    }
    const W normal = random::normal_cosine(d.words);
    const W read = normal * x.source.read_noise_rms_electrons;
    std::optional<double> r, y;
    d.status = store(read, r);
    if (d.status == S::ok) {
      const W hi = W(n) + *r, bp = hi - W(n);
      const W low = (W(n) - (hi - bp)) + (W(*r) - bp);
      d.status = biased_store(hi, low, x.source.gain_electrons_per_adu,
                              x.source.bias_adu, y);
    }
    if (d.status == S::ok) {
      d.poisson_electrons = n;
      d.read_electrons = r;
      d.measured_adu = y;
    }
  }
  return out;
}
} // namespace irred::detector
