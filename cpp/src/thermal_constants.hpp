#pragma once
#include "irred/quantities.hpp"
#include <numbers>
namespace irred::cosmology::detail {
inline constexpr long double thermal_gravitational_constant = 6.67430e-11L;
// Existing fixed-G physical equation/arithmetic sequence, now shared by the
// natural-unit background normalization and pure-H baryon number-density map.
inline long double critical_energy_density_si(double h0_km_s_mpc) {
  const long double c = speed_of_light_m_per_s;
  const long double rate = static_cast<long double>(h0_km_s_mpc) * 1000 /
                           megaparsec_in_metres_wide();
  return 3 * rate * rate * c * c /
         (8 * std::numbers::pi_v<long double> * thermal_gravitational_constant);
}
inline long double critical_mass_density_si(double h0_km_s_mpc) {
  const long double c = speed_of_light_m_per_s;
  return critical_energy_density_si(h0_km_s_mpc) / (c * c);
}
} // namespace irred::cosmology::detail
