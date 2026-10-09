// Independent installed C++ SDK use; synthetic mathematical consumer only.
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <irred/dgp_growth.hpp>
int main() {
  using namespace irred::cosmology;
  auto model = prepare_dgp_growth({.3, 70});
  const std::array<double, 3> a{.1, .5, 1};
  const auto result =
      model.evaluate(a, dgp_e | dgp_h | dgp_omega_m | dgp_mu | dgp_d | dgp_f);
  if (result.status != irred::numerics::Status::ok || result.rows.size() != 3)
    return 1;
  std::cout << std::setprecision(17);
  std::cout << "{\"model\":\"" << dgp_growth_id << "\",\"rows\":[";
  for (size_t i = 0; i < a.size(); ++i) {
    const auto &r = result.rows[i];
    if (!r.e.value || !r.h.value || !r.omega_m.value || !r.mu.value ||
        !r.d.value || !r.f.value)
      return 2;
    if (i)
      std::cout << ',';
    std::cout << "{\"a\":" << a[i] << ",\"E\":" << *r.e.value
              << ",\"H_km_s_Mpc\":" << *r.h.value
              << ",\"Omega_m\":" << *r.omega_m.value
              << ",\"mu\":" << *r.mu.value << ",\"D\":" << *r.d.value
              << ",\"f\":" << *r.f.value << '}';
  }
  std::cout << "],\"interpretation\":\"synthetic dust quasistatic SDK control; "
               "no inference\"}\n";
}
