#pragma once
#include "hydrogen_quantum_density.hpp"
#include "irred/quantities.hpp"
#include <cmath>
namespace irred::cosmology::detail {
struct HydrogenRates {
  long double alpha, beta, alpha_temperature_derivative,
      beta_temperature_derivative;
};
// Original SSS1999 eq3 and p4 detailed-balance relation; F=1.
// B2 is the explicitly declared approximate B1-E21 binding threshold.
inline HydrogenRates hydrogen_rates(long double temperature, long double B1,
                                    long double E21) {
  const long double t = temperature / 10000;
  const long double alpha = 1e-19L * 4.309L * std::pow(t, -.6166L) /
                            (1 + .6703L * std::pow(t, .5300L));
  const long double beta =
      alpha *
      atomic::detail::electron_quantum_density_si(
          boltzmann_constant_joule_per_kelvin * temperature) *
      std::exp(-(B1 - E21) /
               (boltzmann_constant_joule_per_kelvin * temperature));
  const long double logarithmic_alpha =
      (-.6166L - .6703L * .5300L * std::pow(t, .5300L) /
                     (1 + .6703L * std::pow(t, .5300L))) /
      temperature;
  return {alpha, beta, alpha * logarithmic_alpha,
          beta * (logarithmic_alpha + 1.5L / temperature +
                  (B1 - E21) / (boltzmann_constant_joule_per_kelvin *
                                temperature * temperature))};
}
// Exact finite-cell primitive for a linearly interpolated opacity (or its
// nonnegative absolute-error field); theta is measured upward from low z.
inline long double hydrogen_opacity_cell_integral(long double low,
                                                  long double high,
                                                  long double width,
                                                  long double theta) {
  return width * theta * (low + (high - low) * theta / 2);
}
} // namespace irred::cosmology::detail
