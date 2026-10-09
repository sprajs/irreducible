#include "irred/quintessence.hpp"
#include <iostream>
int main() {
  using namespace irred::cosmology;
  auto background=prepare_quintessence({1.2,67.4,.08,.65,.00009});
  const double a[]{1,.8,.5};
  auto b=background.evaluate(a,true);
  if(b.status!=irred::numerics::Status::ok) return 1;
  for(const auto &row:b.rows) {
    if(!row.e.value||!row.w_phi.value||!row.dm_mpc.value) return 2;
    std::cout<<row.scale_factor<<' '<<*row.e.value<<' '<<*row.w_phi.value<<' '<<*row.dm_mpc.value<<'\n';
  }
}
