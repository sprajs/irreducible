#include "irred/thin_lens_geometry.hpp"
#include "irred/effective_neutrino.hpp"
#include "fixtures/thin_lens_geometry_reference.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
int main() {
  using namespace irred::cosmology;using S=irred::numerics::Status;using W=long double;
  std::vector<ThinLensGeometry> owners;
  for(auto s:{CurvedFLRWSpec{70,1,0,0},{70,0,0,0},{67.4,.3,.0001,.5999},
              {67.4,.3,.0001,.7999},{70,10,0,0}})owners.emplace_back(CurvedFLRW(s));
  const auto source=prepare_effective_neutrino_state({67.4,.02237,.12,2.7255,3.046,{{.06,.71611,1}}});
  const auto request=source.observable_request(1059.95,"synthetic supplied drag","matched FD");
  if(!request)return 1;
  owners.emplace_back(prepare_thermal_observables(*request));
  unsigned comparisons=0;
  for(const auto &test:thin_lens_reference::cases) {
    const auto batch=owners[test.model].evaluate({&test.epochs,1},(1u<<lens_geometry_output_count)-1);
    if(batch.status!=S::ok||batch.rows.size()!=1)throw std::runtime_error("reference geometry refusal");
    for(unsigned i=0;i<6;++i) {
      const auto &output=batch.rows[0].outputs[i];
      if(output.status!=S::ok||!output.value)throw std::runtime_error("reference output absent");
      const W absolute=i==3?1e-12L:i==4?1e-10L:1e-8L;
      const W budget=absolute+1e-8L*std::abs(test.values[i]);
      const W reference_error=test.errors[i]+64*std::numeric_limits<W>::epsilon()*std::abs(test.values[i]);
      if(reference_error>.05L*budget)throw std::runtime_error("high-precision reference uncertainty");
      if(std::abs(W(*output.value)-test.values[i])+reference_error>budget)
        throw std::runtime_error("independent high-precision FD/nonflat geometry");
      ++comparisons;
      if(test.has_class) {
        const W class_error=test.class_errors[i]+128*std::numeric_limits<double>::epsilon()*std::abs(test.class_values[i]);
        if(class_error+reference_error>.05L*budget)throw std::runtime_error("CLASS refinement uncertainty");
        if(std::abs(W(*output.value)-test.class_values[i])+class_error>budget)
          throw std::runtime_error("matched refined CLASS geometry");
        ++comparisons;
      }
    }
  }
  std::cout<<"PASS independent lens geometry "<<comparisons<<" comparisons\n";
}
