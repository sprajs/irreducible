// Independently derived observer-frequency quadrature and joint weighted
// moments. Shared SI definitions/interpolation model are ancestry, not
// independent physics. Frozen means1e-280+2e-10rel,
// covariance1e-280+2e-8sqrt(CiiCjj); refinement<=5%. Positive normal references
// cannot be admitted as zero.
#include "irred/photometry_calibration.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::photometry;
namespace {
using W = long double;
unsigned checks = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
W interp(std::span<const double> x, std::span<const double> y, W q) {
  if (q < x.front() || q > x.back())
    return 0;
  size_t j = 1;
  while (j + 1 < x.size() && q > x[j])
    ++j;
  return y[j - 1] +
         (q - x[j - 1]) * (W(y[j]) - y[j - 1]) / (W(x[j]) - x[j - 1]);
}
template <class F> W simpson(F f, W lo, W hi, unsigned n) {
  W sum = f(lo) + f(hi);
  for (unsigned k = 1; k < n; ++k)
    sum += (k % 2 ? 4 : 2) * f(lo + (hi - lo) * k / n);
  return sum * (hi - lo) / (3 * n);
}
std::vector<W> frequency(const CalibrationInput &in,
                         const CalibrationState &state, unsigned n) {
  std::vector<W> v;
  const W unit = 1e-6L, z1 = 1 + W(in.redshift),
          hc = 6.62607015e-34L * 299792458.L,
          scale = 1 / (4 * std::acos(-1.L) * in.luminosity_distance_metre *
                       in.luminosity_distance_metre * z1);
  for (size_t b = 0; b < in.bands.size(); ++b) {
    auto band = in.bands[b];
    auto tr = state.optical_transmission[b];
    W lo = std::max(W(band.observed_wavelength_metre.front()),
                    W(in.spectrum.wavelength_metre.front()) * z1),
      hi = std::min(W(band.observed_wavelength_metre.back()),
                    W(in.spectrum.wavelength_metre.back()) * z1);
    std::vector<W> knots{unit / lo, unit / hi};
    for (double x : band.observed_wavelength_metre)
      if (x > lo && x < hi)
        knots.push_back(unit / x);
    for (double x : in.spectrum.wavelength_metre)
      if (x * z1 > lo && x * z1 < hi)
        knots.push_back(unit / (x * z1));
    std::sort(knots.begin(), knots.end());
    W energy = 0, photons = 0;
    for (size_t j = 1; j < knots.size(); ++j) {
      auto common = [&](W u) {
        W lambda = unit / u;
        return interp(in.spectrum.wavelength_metre,
                      in.spectrum.luminosity_watt_per_metre, lambda / z1) *
               interp(band.observed_wavelength_metre, tr, lambda) * scale;
      };
      energy += simpson([&](W u) { return common(u) * unit / (u * u); },
                        knots[j - 1], knots[j], n);
      photons += simpson(
          [&](W u) { return common(u) * unit * unit / (hc * u * u * u); },
          knots[j - 1], knots[j], n);
    }
    W exposure =
        W(band.collecting_area_square_metre) * band.observer_exposure_second;
    v.push_back(energy * exposure);
    v.push_back(photons * exposure);
  }
  return v;
}
struct Moments {
  std::vector<W> mean, cov;
};
Moments moments(const std::vector<std::vector<W>> &v,
                const CalibrationInput &in) {
  size_t m = v[0].size();
  Moments out{std::vector<W>(m), std::vector<W>(m * m)};
  W total = 0;
  for (auto &s : in.states)
    total += s.relative_mass;
  for (size_t i = 0; i < m; ++i)
    for (size_t s = 0; s < v.size(); ++s)
      out.mean[i] += W(in.states[s].relative_mass) / total * v[s][i];
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < m; ++j)
      for (size_t s = 0; s < v.size(); ++s)
        out.cov[i * m + j] += W(in.states[s].relative_mass) / total *
                              (v[s][i] - out.mean[i]) * (v[s][j] - out.mean[j]);
  return out;
}
void near(W x, W y, W budget, const char *s) {
  need(std::abs(x - y) <= budget, s);
  if (y != 0 && std::abs(y) >= std::numeric_limits<double>::min())
    need(x != 0, "positive normal reference retained");
}
} // namespace
int main() {
  try {
    const double sx[]{1e-6, 2e-6, 5e-6}, sy[]{1000, 2100, 700};
    const double grid[]{1.3e-6, 2.7e-6, 5.9e-6};
    const double a[]{.2, .5, .7}, b[]{.8, .6, .3};
    const CalibrationBand bands[]{{"blue", grid, 2, 10}, {"red", grid, 3, 4}};
    const std::span<const double> t0[]{a, b}, t1[]{b, a}, t2[]{a, b};
    const CalibrationState states[]{
        {"s0", 1, t0}, {"s1", 2, t1}, {"s2", 4, t2}};
    CalibrationInput input{{sx, sy},
                           1e8,
                           .2,
                           bands,
                           states,
                           "synthetic spectrum",
                           "synthetic transmissions",
                           "finite joint masses",
                           "all bands shared states"};
    auto result = evaluate_calibration(input);
    need(result.status == irred::numerics::Status::ok && result.moments &&
             result.axes.size() == 4,
         "joint moments available");
    std::vector<std::vector<W>> fine, coarse;
    for (auto &s : states) {
      fine.push_back(frequency(input, s, 2048));
      coarse.push_back(frequency(input, s, 1024));
    }
    auto ref = moments(fine, input), old = moments(coarse, input);
    for (size_t i = 0; i < 4; ++i) {
      W budget = 1e-280L + 2e-10L * std::abs(ref.mean[i]);
      near(ref.mean[i], old.mean[i], .05L * budget,
           "frequency mean refinement");
      near(result.moments->mean[i], ref.mean[i], budget, "frequency mean");
      for (size_t j = 0; j < 4; ++j) {
        W scale = std::sqrt(ref.cov[i * 4 + i] * ref.cov[j * 4 + j]),
          cbudget = 1e-280L + 2e-8L * scale;
        near(ref.cov[i * 4 + j], old.cov[i * 4 + j], .05L * cbudget,
             "covariance reference refinement");
        near(result.moments->covariance[i * 4 + j], ref.cov[i * 4 + j], cbudget,
             "joint covariance");
      }
    }
    need(ref.cov[2] < 0, "physical anticorrelated bands");
    need(result.axes[0].band_id == "blue" &&
             result.axes[1].output == transmitted_photons &&
             result.attempts.size() == 6,
         "axis and attempt order");
    need(std::abs(result.moments->covariance[5] - result.moments->mean[1]) > 1,
         "calibration signal variance differs from shot noise");
    CalibrationState reordered[]{states[2], states[0], states[1]};
    auto permuted_input = input;
    permuted_input.states = reordered;
    auto permuted = evaluate_calibration(permuted_input);
    need(permuted.moments.has_value(), "state permutation admitted");
    for (size_t i = 0; i < 4; ++i)
      near(permuted.moments->mean[i], result.moments->mean[i],
           1e-280L + 2e-10L * std::abs(ref.mean[i]), "state order invariance");
    CalibrationState scaled[]{states[0], states[1], states[2]};
    for (auto &s : scaled)
      s.relative_mass *= 7;
    permuted_input.states = scaled;
    permuted = evaluate_calibration(permuted_input);
    need(permuted.moments.has_value(), "relative mass scale admitted");
    for (size_t i = 0; i < 16; ++i)
      near(permuted.moments->covariance[i], result.moments->covariance[i],
           1e-280L + 2e-8L * std::sqrt(ref.cov[(i / 4) * 4 + i / 4] *
                                       ref.cov[(i % 4) * 4 + i % 4]),
           "mass normalization scale invariance");
    // Borrowed arrays are synchronous; results own origin/ID strings and
    // moments.
    auto kept = result;
    input.source_origin = "mutated caller";
    need(kept.source_origin == "synthetic spectrum", "owned origins");
    CalibrationPolicy p;
    p.requested_outputs = collected_energy;
    auto energy = evaluate_calibration(input, p);
    need(energy.moments && energy.axes.size() == 2, "omitted photon axis");
    p.maximum_states = 2;
    need(evaluate_calibration(input, p).status ==
             irred::numerics::Status::work_limit,
         "state quota");
    CalibrationState bad[]{states[0], states[1]};
    bad[0].relative_mass = std::numeric_limits<double>::denorm_min();
    input.states = bad;
    need(!evaluate_calibration(input).moments, "subnormal mass refused");
    input.states = {states, 1};
    auto single = evaluate_calibration(input);
    need(single.moments && std::all_of(single.moments->covariance.begin(),
                                       single.moments->covariance.end(),
                                       [](double x) { return x == 0; }),
         "single-state constant witness");
    std::cout << "PASS " << checks
              << " photometry calibration frequency peer controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
