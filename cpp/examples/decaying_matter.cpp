#include "irred/decaying_matter.hpp"
#include <array>
#include <iostream>
int main() {
  using namespace irred::cosmology;
  // Supplied finite epoch, not H0. All densities refer to this same epoch.
  DecayingMatterModel m{.1, 2e-17, .2, .7, .05, .05, .5};
  const std::array scales{.1,.2,.4,.8};
  const auto trajectory=evolve_decaying_matter(m,scales);
  if(trajectory.status!=irred::numerics::Status::ok) return 1;
  std::cout<<decaying_matter_model_id<<"\na E H[s^-1] daughter_fraction q\n";
  for(const auto &r:trajectory.rows)
    std::cout<<r.scale_factor<<' '<<r.expansion_over_initial_hubble<<' '
             <<r.hubble_per_second<<' '<<r.daughter_radiation_fraction<<' '
             <<r.deceleration<<'\n';
}
