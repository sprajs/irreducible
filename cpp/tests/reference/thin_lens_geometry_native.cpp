// Original witness driver, not a production configuration route.
#include "irred/thin_lens_geometry.hpp"
#include "irred/effective_neutrino.hpp"
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>
int main() {
  using namespace irred::cosmology;
  std::vector<std::pair<const char*,ThinLensGeometry>> owners;
  owners.emplace_back("eds",ThinLensGeometry(CurvedFLRW({70,1,0,0})));
  owners.emplace_back("milne",ThinLensGeometry(CurvedFLRW({70,0,0,0})));
  owners.emplace_back("open",ThinLensGeometry(CurvedFLRW({67.4,.3,.0001,.5999})));
  owners.emplace_back("closed",ThinLensGeometry(CurvedFLRW({67.4,.3,.0001,.7999})));
  owners.emplace_back("closed-dust",ThinLensGeometry(CurvedFLRW({70,10,0,0})));
  auto state=prepare_effective_neutrino_state({67.4,.02237,.12,2.7255,3.046,{{.06,.71611,1}}});
  auto request=state.observable_request(1059.95,"synthetic supplied drag","matched FD source");
  if(!request)return 1;
  owners.emplace_back("fd",ThinLensGeometry(prepare_thermal_observables(*request)));
  std::puts("[");bool comma=false;
  for(const auto &[name,owner]:owners) {
    std::vector<LensEpochPair> pairs{{.1,.5},{.5,2},{2,10},{1,std::nextafter(1.,2.)}};
    if(std::string_view(name)=="closed-dust")pairs.push_back({.5,100});
    const auto batch=owner.evaluate(pairs,(1u<<lens_geometry_output_count)-1);
    for(const auto &row:batch.rows) {
      if(comma)std::puts(",");comma=true;
      std::printf("{\"model\":\"%s\",\"lens_redshift\":%.17g,\"source_redshift\":%.17g,\"status\":%u,\"callbacks\":%zu,\"values\":[",name,
          row.epochs.lens_redshift,row.epochs.source_redshift,unsigned(row.status),row.callbacks);
      for(unsigned i=0;i<6;++i) {
        if(i)std::printf(",");if(row.outputs[i].value)std::printf("%.17g",*row.outputs[i].value);else std::printf("null");
      }
      std::printf("],\"errors\":[");for(unsigned i=0;i<6;++i)std::printf("%s%.17g",i?",":"",row.outputs[i].error_estimate);
      std::printf("]}");
    }
  }
  std::puts("\n]");return 0;
}
