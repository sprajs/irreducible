#include <irred/plummer_sphere.hpp>
#include <cmath>
#include <cstdlib>
int main() {
  irred::gravity::PlummerSphere sphere{2e41,3e19,6.67430e-11};
  auto e=irred::gravity::evaluate_plummer_sphere(sphere,3e19);
  auto p=irred::gravity::project_plummer_sphere(sphere,3e19);
  auto l=irred::gravity::plummer_thin_lens(sphere,3e19,20,3e24);
  if(e.status!=irred::numerics::Status::ok || p.status!=irred::numerics::Status::ok || l.status!=irred::numerics::Status::ok) return EXIT_FAILURE;
  if(std::abs(p.projected_mass_kg.value/1e41-1)>1e-14) return EXIT_FAILURE;
  if(!(p.projected_mass_kg.value>e.enclosed_mass_kg.value)) return EXIT_FAILURE;
  if(std::abs(l.convergence.value/(p.surface_density_kg_m2.value/20)-1)>1e-14) return EXIT_FAILURE;
  auto center=irred::gravity::plummer_thin_lens(sphere,0,20,3e24);
  if(center.status!=irred::numerics::Status::ok || center.deflection_radians.value!=0) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
