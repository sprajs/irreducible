#include "irred/photometry_calibration.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
namespace irred::photometry {
namespace {
using S = numerics::Status;
constexpr std::size_t hard_states = 1024, hard_bands = 64, hard_knots = 1048576;
bool add(std::size_t &n, std::size_t x) {
  if (x > std::numeric_limits<std::size_t>::max() - n)
    return false;
  n += x;
  return true;
}
bool mul(std::size_t a, std::size_t b, std::size_t &c) {
  if (b && a > std::numeric_limits<std::size_t>::max() / b)
    return false;
  c = a * b;
  return true;
}
bool cast(long double x, double &y) {
  y = static_cast<double>(x);
  return std::isfinite(x) && std::isfinite(y) && (x == 0 || std::isnormal(y));
}
} // namespace
CalibrationResult evaluate_calibration(const CalibrationInput &in,
                                       CalibrationPolicy p) {
  CalibrationResult out;
  auto fail = [&](S s) {
    out.status = s;
    out.moments.reset();
    return std::move(out);
  };
  const auto ns = in.states.size(), nb = in.bands.size();
  if (std::numeric_limits<long double>::digits < 64)
    return fail(S::conditioning_budget_exceeded);
  const auto mask = p.requested_outputs;
  if (!mask || (mask & ~(collected_energy | transmitted_photons)) ||
      !std::isfinite(p.mean_relative_sensitivity) ||
      !std::isfinite(p.covariance_relative_sensitivity) ||
      p.mean_relative_sensitivity <= 0 ||
      p.covariance_relative_sensitivity <= 0)
    return fail(S::invalid_input);
  if (!ns || !nb)
    return fail(S::invalid_input);
  if (in.spectrum.wavelength_metre.size() !=
      in.spectrum.luminosity_watt_per_metre.size())
    return fail(S::invalid_input);
  if (ns > hard_states || nb > hard_bands || ns > p.maximum_states ||
      nb > p.maximum_bands)
    return fail(S::work_limit);
  std::size_t axes = nb * ((mask & collected_energy ? 1 : 0) +
                           (mask & transmitted_photons ? 1 : 0)),
              attempts = 0, cells = 0, values = 0,
              bytes = sizeof(CalibrationResult),
              knots = in.spectrum.wavelength_metre.size();
  if (!mul(ns, nb, attempts) || !mul(axes, axes, cells) ||
      !mul(ns, axes, values))
    return fail(S::work_limit);
  // Sum declared array lengths before reading samples. Include source each
  // call.
  for (const auto &b : in.bands) {
    std::size_t k = in.spectrum.wavelength_metre.size();
    if (!add(k, b.observed_wavelength_metre.size()) || k > 65536 ||
        !mul(ns, k, k) || !add(knots, k))
      return fail(S::work_limit);
  }
  if (knots > hard_knots || knots > p.maximum_total_knots)
    return fail(S::work_limit);
  auto account = [&](std::size_t count, std::size_t size) {
    std::size_t q = 0;
    return mul(count, size, q) && add(bytes, q);
  };
  if (!account(attempts, sizeof(CalibrationAttempt)) ||
      !account(axes, sizeof(CalibrationAxis) + sizeof(bool) +
                         5 * sizeof(long double) + 2 * sizeof(double)) ||
      !account(cells, 2 * sizeof(long double) + 2 * sizeof(double)) ||
      !account(values, 2 * sizeof(long double)) ||
      !account(ns, sizeof(long double)))
    return fail(S::work_limit);
  for (auto origin : {in.source_origin, in.calibration_origin,
                      in.distribution_origin, in.dependence_origin}) {
    if (origin.empty())
      return fail(S::invalid_input);
    if (origin.size() >= p.maximum_output_bytes)
      return fail(S::work_limit);
    if (!add(bytes, origin.size() + 1))
      return fail(S::work_limit);
  }
  for (const auto &b : in.bands) {
    if (b.id.size() >= p.maximum_output_bytes)
      return fail(S::work_limit);
    std::size_t size = 0;
    if (!mul(ns + 2, b.id.size() + 1, size) || !add(bytes, size))
      return fail(S::work_limit);
  }
  for (const auto &s : in.states) {
    if (s.id.size() >= p.maximum_output_bytes)
      return fail(S::work_limit);
    std::size_t size = 0;
    if (!mul(nb, s.id.size() + 1, size) || !add(bytes, size))
      return fail(S::work_limit);
  }
  if (bytes > p.maximum_output_bytes)
    return fail(S::work_limit);
  try {
    // Quadratic duplicate checking is bounded and needs no additional heap.
    for (std::size_t i = 0; i < nb; ++i) {
      if (in.bands[i].id.empty())
        return fail(S::invalid_input);
      for (std::size_t j = 0; j < i; ++j)
        if (in.bands[i].id == in.bands[j].id)
          return fail(S::invalid_input);
    }
    long double mass = 0;
    for (std::size_t s = 0; s < ns; ++s) {
      const auto &st = in.states[s];
      if (st.id.empty() || !std::isnormal(st.relative_mass) ||
          st.relative_mass <= 0 || st.optical_transmission.size() != nb)
        return fail(S::invalid_input);
      for (std::size_t t = 0; t < s; ++t)
        if (st.id == in.states[t].id)
          return fail(S::invalid_input);
      for (std::size_t b = 0; b < nb; ++b)
        if (st.optical_transmission[b].size() !=
            in.bands[b].observed_wavelength_metre.size())
          return fail(S::invalid_input);
      mass += static_cast<long double>(st.relative_mass);
    }
    if (!std::isfinite(mass) || mass <= 0)
      return fail(S::overflow);
    std::vector<long double> weights(ns), f(values), err(values), mean(axes),
        me(axes), cov(cells), ce(cells);
    std::vector<bool> constant(axes, true);
    out.axes.reserve(axes);
    out.attempts.reserve(attempts);
    out.source_origin = in.source_origin;
    out.calibration_origin = in.calibration_origin;
    out.distribution_origin = in.distribution_origin;
    out.dependence_origin = in.dependence_origin;
    bool source_zero =
        !in.spectrum.luminosity_watt_per_metre.empty() &&
        std::all_of(in.spectrum.luminosity_watt_per_metre.begin(),
                    in.spectrum.luminosity_watt_per_metre.end(),
                    [](double x) { return x == 0; });
    std::size_t a = 0;
    for (std::size_t b = 0; b < nb; ++b) {
      bool fixed = true;
      for (std::size_t s = 1; s < ns; ++s)
        if (!std::equal(in.states[0].optical_transmission[b].begin(),
                        in.states[0].optical_transmission[b].end(),
                        in.states[s].optical_transmission[b].begin()))
          fixed = false;
      fixed = fixed || source_zero ||
              in.bands[b].collecting_area_square_metre == 0 ||
              in.bands[b].observer_exposure_second == 0;
      for (auto output : {collected_energy, transmitted_photons})
        if (mask & output) {
          out.axes.push_back({std::string(in.bands[b].id), output});
          constant[a++] = fixed;
        }
    }
    S failure = S::ok;
    const auto eps = std::numeric_limits<long double>::epsilon();
    for (std::size_t s = 0; s < ns; ++s) {
      weights[s] = static_cast<long double>(in.states[s].relative_mass) / mass;
      a = 0;
      for (std::size_t b = 0; b < nb; ++b) {
        const auto &band = in.bands[b];
        SampledInput input{in.spectrum,
                           {band.observed_wavelength_metre,
                            in.states[s].optical_transmission[b]},
                           in.luminosity_distance_metre,
                           in.redshift,
                           band.collecting_area_square_metre,
                           band.observer_exposure_second};
        auto result = evaluate_sampled(input, {mask, 65536, 65536});
        out.attempts.push_back(
            {std::string(in.states[s].id), std::string(band.id), result});
        if (result.admission_status != S::ok)
          failure = result.admission_status;
        for (auto output : {collected_energy, transmitted_photons})
          if (mask & output) {
            const auto &o = output == collected_energy
                                ? result.energy_joule
                                : result.expected_photons;
            if (o.availability != Availability::available || !o.value ||
                o.numerical_status != S::ok) {
              failure = o.numerical_status == S::ok ? result.admission_status
                                                    : o.numerical_status;
            } else {
              f[s * axes + a] = *o.value;
              err[s * axes + a] = 3e-12L * std::abs(*o.value);
              mean[a] += weights[s] * f[s * axes + a];
              me[a] += weights[s] * err[s * axes + a];
            }
            ++a;
          }
      }
    }
    if (failure != S::ok)
      return fail(failure);
    for (a = 0; a < axes; ++a) {
      me[a] += (ns + 4) * eps * std::abs(mean[a]);
      if (mean[a] > 0 && me[a] > p.mean_relative_sensitivity * mean[a])
        return fail(S::conditioning_budget_exceeded);
    }
    for (std::size_t i = 0; i < axes; ++i)
      for (std::size_t j = 0; j <= i; ++j) {
        long double v = 0, e = 0, absolute_products = 0;
        if (!constant[i] && !constant[j])
          for (std::size_t s = 0; s < ns; ++s) {
            const auto di = f[s * axes + i] - mean[i],
                       dj = f[s * axes + j] - mean[j],
                       ei = err[s * axes + i] + me[i],
                       ej = err[s * axes + j] + me[j];
            v += weights[s] * di * dj;
            absolute_products += weights[s] * std::abs(di * dj);
            e += weights[s] * (std::abs(di) * ej + std::abs(dj) * ei + ei * ej);
          }
        e += (ns + 8) * eps * absolute_products;
        cov[i * axes + j] = cov[j * axes + i] = v;
        ce[i * axes + j] = ce[j * axes + i] = e;
      }
    for (std::size_t i = 0; i < axes; ++i) {
      if (!constant[i] && !(cov[i * axes + i] > 0))
        return fail(S::conditioning_budget_exceeded);
      for (std::size_t j = 0; j < axes; ++j) {
        auto scale =
            std::sqrt(cov[i * axes + i]) * std::sqrt(cov[j * axes + j]);
        if (ce[i * axes + j] > p.covariance_relative_sensitivity * scale)
          return fail(S::conditioning_budget_exceeded);
      }
    }
    CalibrationMoments m;
    m.mean.resize(axes);
    m.mean_empirical_sensitivity.resize(axes);
    m.covariance.resize(cells);
    m.covariance_empirical_sensitivity.resize(cells);
    for (a = 0; a < axes; ++a) {
      if (!cast(mean[a], m.mean[a]))
        return fail(S::overflow);
      double e = 0;
      if (!cast(me[a], e))
        return fail(S::overflow);
      m.mean_empirical_sensitivity[a] = e + std::abs(m.mean[a] - mean[a]);
    }
    for (a = 0; a < cells; ++a) {
      if (!cast(cov[a], m.covariance[a]))
        return fail(S::overflow);
      double e = 0;
      if (!cast(ce[a], e))
        return fail(S::overflow);
      m.covariance_empirical_sensitivity[a] =
          e + std::abs(m.covariance[a] - cov[a]);
    }
    out.moments = std::move(m);
    out.status = S::ok;
    return out;
  } catch (const std::bad_alloc &) {
    return fail(S::work_limit);
  } catch (const std::length_error &) {
    return fail(S::work_limit);
  }
}
} // namespace irred::photometry
