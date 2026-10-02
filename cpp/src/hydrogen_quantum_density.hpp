#pragma once
#include "irred/hydrogen_equilibrium.hpp"
#include "irred/quantities.hpp"
#include <cmath>
#include <numbers>
namespace irred::atomic::detail {
// Shared dilute electron translational quantum density. The scalar expression,
// fixed mass and qbase*sqrt(qbase) arithmetic preserve the original Saha
// source.
inline long double electron_quantum_density_si(long double kT) {
  const long double qbase =
      2 * std::numbers::pi_v<long double> * hydrogen_electron_mass_kg * kT /
      (planck_constant_joule_second * planck_constant_joule_second);
  return qbase * std::sqrt(qbase);
}
} // namespace irred::atomic::detail
