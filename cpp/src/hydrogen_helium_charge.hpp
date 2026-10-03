#pragma once
#include "irred/quantities.hpp"
namespace irred::cosmology::detail {
// Shared two-stage charge and fixed Thomson asset. Preserve the existing
// multiplication/addition sequence; this does not add a HeIII population.
inline constexpr long double hhe_thomson_cross_section_si = 6.6524587051e-29L;
inline long double hhe_electron_density(long double nH, long double p,
                                        long double nHe, long double q) noexcept {
  return nH * p + nHe * q;
}
inline long double hhe_opacity_coefficient(long double H,
                                          long double one_plus_z) noexcept {
  return speed_of_light_m_per_s * hhe_thomson_cross_section_si /
         (H * one_plus_z);
}
} // namespace irred::cosmology::detail
