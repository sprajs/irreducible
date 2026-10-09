#include "irred/chaplygin.hpp"
#include <array>
#include <iostream>
int main(){
 using namespace irred::cosmology;
 // Ordinary matter is separate from the unified fluid; Omega_cg=.85 is derived.
 ChaplyginModel model{70,.1,.05,.75,.8};
 const std::array scales{1.,.7,.3,.1};
 const auto trajectory=evolve_chaplygin(model,scales,chaplygin_distances);
 if(trajectory.status!=irred::numerics::Status::ok)return 1;
 std::cout<<chaplygin_model_id<<"\na H[km/s/Mpc] w_cg formal_dp/de D_L[Mpc]\n";
 for(const auto &r:trajectory.rows)
  std::cout<<r.scale_factor<<' '<<r.hubble_km_s_mpc<<' '<<r.fluid_equation_of_state<<' '
           <<r.fluid_barotropic_slope<<' '<<r.distances->luminosity_mpc<<'\n';
}
