#pragma once
#include <cmath>
namespace irred::cosmology::detail {
// Shared radiation-free flat matter+Lambda equation. Public late Expansion
// keeps its own existing admission. Consumers must qualify their source/domain.
inline long double lcdm_lambda(double omega_m) {
  return 1 - (long double)omega_m;
}
inline long double lcdm_expansion(double omega_m, long double z) {
  const auto u = 1 + z;
  return std::sqrt(omega_m * u * u * u + lcdm_lambda(omega_m));
}
inline long double lcdm_scaled_expansion(double omega_m, long double a) {
  return omega_m + lcdm_lambda(omega_m) * a * a * a;
}
} // namespace irred::cosmology::detail
