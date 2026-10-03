#include <irred/helium_escape_history.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <utility>
int main() {
  using namespace irred::cosmology;
  using S = irred::numerics::Status;
  HeliumEscapeDriverKnot rows[]{ {2700, 1e-13, 3e9, 6000, 1e-5},
                                {1900, 1e-13, 3e9, 6000, 1e-5} };
  const long double f = double(.08), p = 1 - static_cast<long double>(rows[0].hydrogen_neutral_fraction),
      s = 4 * 2.414194e21L * 6000 * std::sqrt(6000.L) / 3e9L * std::exp(-285325.L / 6000),
      q = 2 * s / (p + s + std::sqrt((p + s) * (p + s) + 4 * f * s));
  HeliumEscapeHistoryRequest source{rows, .08, static_cast<double>(q),
      "installed SDK / synthetic constant bath", HeliumEscapeDriverRole::synthetic_prescribed_bath};
  auto owner = prepare_helium_escape_history(source);
  if (owner.status() != S::ok || !owner.source()) throw std::runtime_error("installed prepare");
  rows[0].radiation_temperature_kelvin = 1;
  auto copied = owner; auto moved = std::move(copied);
  const double z[]{2700, 2400, 1900}; const auto result = moved.evaluate(z, 7);
  std::cout << "installed native_model=" << helium_escape_history_model_id
            << " native_method=" << helium_escape_history_method_id
            << " work=" << moved.work().total() << " driver_queries=" << result.driver_evaluations << '\n';
  if (result.status != S::ok || result.rows.size() != 3 || copied.source() ||
      moved.source()->knots[0].radiation_temperature_kelvin != 6000)
    throw std::runtime_error("installed ownership/query");
  for (const auto &r : result.rows)
    if (!r.helium_singly_ionized_fraction.value || !r.electron_number_density_per_cubic_metre.value ||
        !r.thomson_opacity_per_redshift.value) throw std::runtime_error("installed requested group");
}
