#pragma once
#include "hydrogen_quantum_density.hpp"
#include "hydrogen_helium_charge.hpp"
#include <cmath>
#include <numbers>
namespace irred::cosmology::detail {
// Independently authored equations from the frozen primary reconciliation.
// Finite scalar facts retain CLASS0ceb/HyRec2020 helium source identity. No
// external implementation/table has been copied. All constants alpha/me=1.
inline long double helium_escape_quantum_density(long double T) {
  const long double k = boltzmann_constant_joule_per_kelvin;
  return 2.414194e21L * atomic::detail::electron_quantum_density_si(k * T) /
         atomic::detail::electron_quantum_density_si(k);
}
inline long double helium_escape_depth_per_neutral(long double H, long double nH,
                                                   long double f) {
  return 4.277e-14L * nH * f / H;
}
inline long double helium_escape_optical_depth(long double H, long double nH,
                                              long double f, long double q) {
  return helium_escape_depth_per_neutral(H, nH, f) * (1 - q);
}
struct HeliumEscapeRate {
  long double tau = 0, continuum_inverse_seconds = 0, width_per_second = 0,
      continuum_depth = 0, enhancement = 0, effective_escape = 0,
      downward_per_second = 0, saha = 0, electron_density = 0,
      value = 0, fraction_derivative = 0, absolute_rate_scale = 0;
};
// Caller admits source/profile and q. Derivative is in ell=ln(a), per He.
// No H derivative term belongs here; ne uses the same shared charge owner.
inline HeliumEscapeRate helium_escape_rate(long double H, long double nH,
                                           long double T, long double xHI,
                                           long double f, long double q) {
  const long double nHe = nH * f,
      ne = hhe_electron_density(nH, 1 - xHI, nHe, q),
      xe = ne / nH, fx = nHe / nH,
      s0 = 4 * helium_escape_quantum_density(T) / nH,
      s = s0 * std::exp(-285325.L / T),
      ys = std::exp(46090.L / T) / s0,
      yp = 3 * std::exp(39101.L / T) / s0,
      inverse = 9.15776e28L * H / (nH * xHI),
      width = 1.976e6L / -std::expm1(-6989.L / T) +
          6.03e6L / std::expm1(19754.L / T) +
          1.06e8L / std::expm1(21539.L / T) +
          2.18e6L / std::expm1(28496.L / T) +
          3.37e7L / std::expm1(29224.L / T) +
          1.04e6L / std::expm1(32414.L / T) +
          1.51e7L / std::expm1(32781.L / T),
      depth_per_neutral = helium_escape_depth_per_neutral(H, nH, f),
      tau = helium_escape_optical_depth(H, nH, f, q),
      k = width / (4 * std::numbers::pi_v<long double> *
                       std::numbers::pi_v<long double> * inverse),
      tc = k * tau, pi2 = std::numbers::pi_v<long double> *
                          std::numbers::pi_v<long double>,
      root = std::sqrt(1 + pi2 * tc),
      enhancement = root + 7.74L * tc / (1 + 70 * tc),
      enhancement_tau = k * (pi2 / (2 * root) +
                                  7.74L / ((1 + 70 * tc) * (1 + 70 * tc))),
      r = 1.023e-7L, survival = std::exp(-r * tau),
      B = -std::expm1(-r * tau), Bt = r * survival,
      phi = 6.14e13L / inverse, D = -std::expm1(-phi),
      R = .964525L * std::exp(2947.L / T),
      numerator = enhancement * (survival + B * D) + B * R,
      numerator_tau = enhancement_tau * (survival + B * D) +
                      Bt * (R - enhancement * std::exp(-phi)),
      escape = numerator / tau,
      escape_q = -depth_per_neutral *
                 (numerator_tau * tau - numerator) / (tau * tau),
      down = 50.94L * ys + 1.7989e9L * yp * escape,
      down_q = 1.7989e9L * yp * escape_q,
      residual = (1 - q) * s - q * xe,
      residual_q = -s - xe - q * fx;
  return {tau, inverse, width, tc, enhancement, escape, down, s, ne,
          down * residual / H,
          (down_q * residual + down * residual_q) / H,
          down * ((1 - q) * s + q * xe) / H};
}
inline long double helium_escape_heiii_log_activity(long double T,
                                                    long double ne) {
  return std::log(helium_escape_quantum_density(T)) - 631462.7L / T -
         std::log(ne);
}
} // namespace irred::cosmology::detail
