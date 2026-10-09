#include <irred/nfw_halo.hpp>
#include <cmath>
#include <cstdlib>
int main() {
  using namespace irred::lensing;
  NFWHalo halo{1e15,0.2};
  auto mass=nfw_enclosed_mass(halo,0.2);
  auto projection=nfw_projection(halo,0.2);
  auto lens=nfw_lens(halo,0.2,3e15,1000);
  if(mass.status!=irred::numerics::Status::ok || projection.status!=irred::numerics::Status::ok || lens.status!=irred::numerics::Status::ok) return EXIT_FAILURE;
  if(std::abs(lens.convergence.value-projection.surface_density.value/3e15)>1e-14) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
