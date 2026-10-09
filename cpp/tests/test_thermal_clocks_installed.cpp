// This consumer is also compiled separately against the installed archive.
#include <irred/effective_neutrino.hpp>
#include <irred/thermal_clocks.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
int main(int argc,char **) {
  using namespace irred::cosmology;
  const auto source=prepare_effective_neutrino_state({67.4,.02237,.12,2.7255,3.046,{{.06,.71611,1}}});
  const auto clocks=prepare_thermal_clocks(source.background(),{},source.diagnostics().stress_source_estimate);
  const double a[]{.001,.01,.1,.5,1,std::nextafter(1.,0.)};
  const auto result=clocks.evaluate(a,(1u<<thermal_clock_output_count)-1);
  if(source.status()!=irred::numerics::Status::ok||clocks.status()!=irred::numerics::Status::ok||result.rows.size()!=6)return 1;
  for(const auto &row:result.rows) {
    for(unsigned i=0;i<row.outputs.size();++i) {
      const auto &out=row.outputs[i];
      if(out.status!=irred::numerics::Status::ok||!out.value) {
        std::cerr<<"clock refusal a="<<row.scale_factor<<" output="<<i
                 <<" status="<<static_cast<unsigned>(out.status)<<" callbacks="<<row.callbacks
                 <<" outer="<<row.outer_callbacks<<" momentum="<<row.momentum_callbacks<<'\n';
        return 2;
      }
    }
    if(argc>1) {
      std::cout<<std::setprecision(17)<<row.scale_factor;
      for(const auto &out:row.outputs)std::cout<<' '<<*out.value<<' '<<out.error_estimate;
      std::cout<<' '<<row.callbacks<<'\n';
    }
  }
  if(argc==1)std::cout<<"PASS installed thermal clock consumer\n";
}
