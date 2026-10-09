#include "irred/bianchi_i.hpp"
#include <array>
#include <iostream>
int main() {
 using namespace irred::cosmology;
 // sx=.3, sy=-.1, sz=-.2 -> Omega_shear=.14/6.
 // All supplied fluid fractions refer to the same a=1 observer anchor.
 BianchiIModel model{2.2e-18,.3,.1,1-.3-.1-.14/6,.3,-.1};
 const std::array<BianchiIQuery,4> queries{{{1,{1,0,0}},
  {.5,{1,0,0}},{.5,{0,1,0}},{.2,{.6,.8,0}}}};
 const auto trajectory=evolve_bianchi_i(model,queries);
 if(trajectory.status!=irred::numerics::Status::ok)return 1;
 std::cout<<bianchi_i_model_id<<"\na E directional_z shear_fraction\n";
 for(const auto &r:trajectory.rows)
  std::cout<<r.query.scale_factor<<' '<<r.expansion_over_anchor_hubble<<' '
           <<r.directional_redshift<<' '<<r.shear_fraction<<'\n';
}
