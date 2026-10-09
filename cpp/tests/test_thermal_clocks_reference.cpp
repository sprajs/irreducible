#include "irred/effective_neutrino.hpp"
#include "irred/thermal_clocks.hpp"
#include "fixtures/thermal_clocks_reference.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
int main() {
  try {
    using namespace irred::cosmology;using S=irred::numerics::Status;
    const auto source=prepare_effective_neutrino_state({67.4,.02237,.12,2.7255,3.046,{{.06,.71611,1}}});
    const auto clocks=prepare_thermal_clocks(source.background(),{},source.diagnostics().stress_source_estimate);
    const auto batch=clocks.evaluate(thermal_clock_reference::factors,(1u<<thermal_clock_output_count)-1);
    if(source.status()!=S::ok||clocks.status()!=S::ok||batch.rows.size()!=6)throw std::runtime_error("reference source admission");
    unsigned checks=0;
    for(unsigned i=0;i<batch.rows.size();++i)for(unsigned j=0;j<thermal_clock_output_count;++j) {
      const auto &out=batch.rows[i].outputs[j];
      if(out.status!=S::ok||!out.value)throw std::runtime_error("reference clock refused");
      const auto reference=thermal_clock_reference::high_precision[i][j];
      const auto budget=1e-8L+2e-10L*std::abs(reference);
      if(std::abs(static_cast<long double>(*out.value)-reference)>budget)throw std::runtime_error("independent full-support FD clock");
      ++checks;
      if(i<thermal_clock_reference::matched_class.size()) {
        if(std::abs(static_cast<long double>(*out.value)-thermal_clock_reference::matched_class[i][j])>budget)
          throw std::runtime_error("independent matched CLASS clock");
        ++checks;
      }
    }
    std::cout<<"PASS independent thermal clock references "<<checks<<" checks\n";
  } catch(const std::exception &e) {std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
