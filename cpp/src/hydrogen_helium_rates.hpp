#pragma once
#include "hydrogen_rates.hpp"
#include "hydrogen_helium_charge.hpp"
#include "irred/hydrogen_helium_equilibrium.hpp"
#include <array>
#include <numbers>
namespace irred::cosmology::detail {
using HHeWide = long double;
inline constexpr HHeWide hhe_lya = 121.5682e-9L,
    hhe_hydrogen_two_photon = 8.22458L, hhe_helium_two_photon = 51.3L,
    hhe_helium_2s_inverse_metres = 1.66277434e7L,
    hhe_helium_2p_inverse_metres = 1.71134891e7L,
    hhe_thomson = hhe_thomson_cross_section_si;
struct HHeRates { HHeWide alpha, beta, photo, alpha_T, beta_T, photo_T; };
inline HHeRates hhe_rates(HHeWide temperature, bool helium) {
  const HHeWide k = boltzmann_constant_joule_per_kelvin,
      chi = electron_volt_joule *
            (helium ? atomic::helium_first_ionization_energy_ev
                    : atomic::hydrogen_ionization_energy_ev),
      excitation = planck_constant_joule_second * speed_of_light_m_per_s *
                   (helium ? hhe_helium_2s_inverse_metres : 1 / hhe_lya);
  HHeWide alpha, alpha_T, beta, beta_T;
  if (!helium) {
    const auto h = hydrogen_rates(temperature, chi, excitation);
    alpha = h.alpha; alpha_T = h.alpha_temperature_derivative;
    beta = h.beta; beta_T = h.beta_temperature_derivative;
  } else {
    const HHeWide s0 = std::sqrt(temperature / std::pow(10.L, .477121L)),
        s1 = std::sqrt(temperature / std::pow(10.L, 5.114L));
    alpha = std::pow(10.L, -16.744L) /
            (s0 * std::pow(1 + s0, 1 - .711L) *
             std::pow(1 + s1, 1 + .711L));
    alpha_T = -alpha * (1 + (1 - .711L) * s0 / (1 + s0) +
                        (1 + .711L) * s1 / (1 + s1)) / (2 * temperature);
    beta = 4 * alpha * atomic::detail::electron_quantum_density_si(k * temperature) *
           std::exp(-(chi - excitation) / (k * temperature));
    beta_T = beta * (alpha_T / alpha + 1.5L / temperature +
                    (chi - excitation) / (k * temperature * temperature));
  }
  // Same detailed-balance equation, combining the exponentials before they
  // underflow; the two-stage stationary law uses exactly the shared chi.
  const HHeWide photo = (helium ? 4 : 1) * alpha *
      atomic::detail::electron_quantum_density_si(k * temperature) *
      std::exp(-chi / (k * temperature));
  return {alpha, beta, photo, alpha_T, beta_T,
          photo * (alpha_T / alpha + 1.5L / temperature +
                   chi / (k * temperature * temperature))};
}
struct HHeRhs {
  std::array<HHeWide, 2> value{}, temperature_derivative{};
  std::array<std::array<HHeWide, 2>, 2> fraction_derivative{};
};
inline HHeRhs hhe_rhs(HHeWide hydrogen, HHeWide helium, HHeWide T,
                     HHeWide nH, HHeWide nHe, HHeWide H, HHeWide u) {
  const HHeWide x[]{hydrogen, helium}, n[]{nH, nHe},
      ne = hhe_electron_density(nH, hydrogen, nHe, helium),
      delta = planck_constant_joule_second * speed_of_light_m_per_s *
              (hhe_helium_2p_inverse_metres - hhe_helium_2s_inverse_metres),
      k = boltzmann_constant_joule_per_kelvin;
  HHeRhs out;
  for (unsigned i = 0; i < 2; ++i) {
    const auto a = hhe_rates(T, i == 1);
    const HHeWide wavelength = i ? 1 / hhe_helium_2p_inverse_metres : hhe_lya,
        K = wavelength * wavelength * wavelength /
            (8 * std::numbers::pi_v<HHeWide> * H),
        L = i ? hhe_helium_two_photon : hhe_hydrogen_two_photon,
        w = i ? std::exp(-delta / (k * T)) : 1,
        wT = i ? w * delta / (k * T * T) : 0,
        v = K * n[i] * (1 - x[i]), numerator = w + L * v,
        denominator = w + (L + a.beta) * v,
        C = numerator / denominator,
        Cx = w * a.beta * K * n[i] / (denominator * denominator),
        CT = (wT * a.beta * v - numerator * a.beta_T * v) /
             (denominator * denominator),
        r = a.alpha * ne * x[i] - a.photo * (1 - x[i]);
    out.value[i] = C * r / (H * u);
    out.temperature_derivative[i] =
        (CT * r + C * (a.alpha_T * ne * x[i] -
                      a.photo_T * (1 - x[i]))) / (H * u);
    for (unsigned j = 0; j < 2; ++j)
      out.fraction_derivative[i][j] =
          (C * (a.alpha * x[i] * n[j] +
                (i == j ? a.alpha * ne + a.photo : 0)) +
           (i == j ? Cx * r : 0)) / (H * u);
  }
  return out;
}
} // namespace irred::cosmology::detail
