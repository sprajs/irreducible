// Independently compiled against installed public headers/archive.
#include <irred/thin_lens_geometry.hpp>
#include <irred/plummer_sphere.hpp>
#include <irred/quantities.hpp>
#include <cmath>
#include <iostream>
int main() {
  using namespace irred::cosmology;using S=irred::numerics::Status;
  ThinLensGeometry geometry(CurvedFLRW({70,.3,0,.7}));
  const LensEpochPair pair{.5,2};
  const unsigned mask=lens_geometry_mask(LensGeometryOutput::lens_angular_mpc)|
      lens_geometry_mask(LensGeometryOutput::critical_surface_density_kg_m2);
  const auto batch=geometry.evaluate({&pair,1},mask);
  if(batch.status!=S::ok)return 1;
  const double dl=*batch.rows[0].outputs[0].value*irred::megaparsec_in_metres_wide();
  const double sigma=*batch.rows[0].outputs[4].value;
  const irred::gravity::PlummerSphere profile{2e41,1e20,6.67430e-11};
  const auto projected=irred::gravity::project_plummer_sphere(profile,2e20);
  const auto lens=irred::gravity::plummer_thin_lens(profile,2e20,sigma,dl);
  if(projected.status!=S::ok||lens.status!=S::ok||sigma<=0||
      std::abs(lens.convergence.value-projected.surface_density_kg_m2.value/sigma)>1e-12)return 1;
  // Geometry diagnostics are deterministic numerical estimates, not profile
  // uncertainties. Propagate Sigma's positive interval through kappa explicitly.
  const double e=batch.rows[0].outputs[4].error_estimate;
  if(!(sigma>e))return 1;
  const double lower=projected.surface_density_kg_m2.value/(sigma+e);
  const double upper=projected.surface_density_kg_m2.value/(sigma-e);
  if(lens.convergence.value<lower||lens.convergence.value>upper)return 1;
  std::cout<<"PASS installed physical-SI geometry to Plummer consumer\n";
}
