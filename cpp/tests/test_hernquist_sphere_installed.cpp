#include <irred/hernquist_sphere.hpp>
#include <cmath>
#include <cstdlib>
int main() {
  // Synthetic SI parameters and explicitly supplied G central value.
  irred::gravity::HernquistSphere sphere{2e41,3e19,6.67430e-11};
  auto e=irred::gravity::evaluate_hernquist_sphere(sphere,3e19);
  if(e.status!=irred::numerics::Status::ok) return EXIT_FAILURE;
  if(std::abs(e.enclosed_mass_kg.value/5e40-1)>1e-14) return EXIT_FAILURE;
  if(std::abs(e.circular_speed_squared_m2_s2.value/(-3e19*e.inward_radial_acceleration_m_s2.value)-1)>1e-14) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
