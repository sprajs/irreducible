#include "irred/curved_flrw.hpp"
#include <stdexcept>
#include <iostream>
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
#define CHECK(expression) require((expression),#expression)
int main(){
  irred::cosmology::CurvedFLRW model({67,.3,.0001,.6});
  const double redshifts[]={.1,.5,1,2};
  auto prediction=model.evaluate(redshifts);
  CHECK(prediction.status==irred::numerics::Status::ok);
  std::cout<<irred::cosmology::curved_flrw_id<<'\n';
  for(const auto& row:prediction.rows)std::cout<<row.redshift<<' '<<row.h_km_s_mpc<<' '<<row.transverse_mpc<<' '<<row.luminosity_mpc<<'\n';
}
