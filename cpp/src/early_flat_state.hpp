#pragma once
#include "irred/sound_horizon.hpp"
#include <algorithm>
#include <cmath>
namespace irred::cosmology::detail {
struct EarlyFlatState {
  long double radiation, matter, lambda, baryon_loading;
};
inline numerics::Status prepare_early_flat(const EarlyFlatModel &m,
                                           EarlyFlatState &s) {
  const double values[]{m.h0_km_s_mpc, m.omega_m, m.omega_r, m.omega_b,
                        m.omega_gamma};
  for (double x : values)
    if (!std::isfinite(x))
      return numerics::Status::nonfinite_input;
  if (!(m.h0_km_s_mpc > 0) || m.omega_m < 0 || !(m.omega_r > 0) ||
      m.omega_b < 0 || m.omega_b > m.omega_m || !(m.omega_gamma > 0) ||
      m.omega_gamma > m.omega_r)
    return numerics::Status::outside_domain;
  // Sterbenz-safe closure rejects positive excess even if the sum rounds to 1.
  const long double largest = std::max(m.omega_m, m.omega_r);
  const long double smaller = std::min(m.omega_m, m.omega_r);
  const long double remaining = 1.L - largest;
  if (largest > 1 || smaller > remaining)
    return numerics::Status::outside_domain;
  s = {(long double)m.omega_r, (long double)m.omega_m, remaining - smaller,
       (3.L * m.omega_b) / (4.L * m.omega_gamma)};
  return numerics::Status::ok;
}
inline long double early_flat_polynomial(const EarlyFlatState &s,
                                         long double a) {
  const long double a2 = a * a;
  return s.radiation + s.matter * a + s.lambda * a2 * a2;
}
} // namespace irred::cosmology::detail
