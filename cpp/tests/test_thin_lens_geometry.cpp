#include "irred/thin_lens_geometry.hpp"
#include "irred/effective_neutrino.hpp"
#include "irred/quantities.hpp"
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
using namespace irred::cosmology;
using S=irred::numerics::Status;
using W=long double;
constexpr unsigned all=(1u<<lens_geometry_output_count)-1;
unsigned checks=0;
void need(bool x,const char *text) { ++checks;if(!x)throw std::runtime_error(text); }
void compare(const LensGeometryRow &row,W dl,W ds,W dls) {
  const W c=irred::speed_of_light_m_per_s,mpc=irred::megaparsec_in_metres_wide();
  const W reference[]{dl,ds,dls,dls/ds,c*c/(4*std::numbers::pi_v<W>*6.67430e-11L)*ds/(dl*dls)/mpc,
                       (1+W(row.epochs.lens_redshift))*dl*ds/dls};
  for(unsigned i=0;i<6;++i) {
    const auto &v=row.outputs[i];
    if(v.status!=S::ok)std::cerr<<"refused output="<<i<<" status="<<unsigned(v.status)<<" work="<<row.callbacks<<'\n';
    need(v.status==S::ok&&v.value.has_value(),"analytic output admitted");
    const W absolute=i==3?1e-12L:i==4?1e-10L:1e-8L;
    const W budget=absolute+1e-8L*std::abs(reference[i]);
    need(std::abs(W(*v.value)-reference[i])<=budget,"independent analytic geometry");
    need(v.error_estimate>=0&&v.error_estimate<=budget,"reported admission");
  }
}
int main() {
  const LensEpochPair pairs[]{{.1,.5},{.5,2},{2,10},{1,std::nextafter(1.,2.)}};
  for(bool milne:{false,true}) {
    ThinLensGeometry owner(CurvedFLRW({70,milne?0.:1.,0,0}));
    need(owner.status()==S::ok&&owner.curved_provider()&&!owner.thermal_provider(),"typed curved owner");
    const auto result=owner.evaluate(pairs,all);need(result.status==S::ok&&result.rows.size()==4,"ordered analytic pairs");
    const W dh=299792.458L/70;
    for(const auto &row:result.rows) {
      const W l=row.epochs.lens_redshift,s=row.epochs.source_redshift;
      // Analytic direct interval uses stable ratio/expm1 for close epochs.
      const W xl=milne?std::log1p(l):2*(1-1/std::sqrt(1+l));
      const W xs=milne?std::log1p(s):2*(1-1/std::sqrt(1+s));
      const W interval=milne?std::log1p((s-l)/(1+l)):
          -2*std::expm1(-.5L*std::log1p((s-l)/(1+l)))/std::sqrt(1+l);
      const W dl=dh*(milne?std::sinh(xl):xl)/(1+l);
      const W ds=dh*(milne?std::sinh(xs):xs)/(1+s);
      const W dls=dh*(milne?std::sinh(interval):interval)/(1+s);
      compare(row,dl,ds,dls);
      need(row.callbacks==row.outer_callbacks+row.momentum_callbacks&&row.momentum_callbacks==0,"curved receipt");
    }
  }
  // Exact closed-dust radial integral; source past pi/2 but before pi.
  const LensEpochPair closed_pair{.5,100};
  const auto closed=ThinLensGeometry(CurvedFLRW({70,10,0,0})).evaluate({&closed_pair,1},all);
  need(closed.status==S::ok,"permitted negative cosine branch");
  const W q=3,dh=299792.458L/70;
  const auto phase=[](W z) { return 2*(std::atan(std::sqrt(1+10*z)/3)-std::atan(1.L/3)); };
  const W pl=phase(closed_pair.lens_redshift),ps=phase(closed_pair.source_redshift);
  need(ps>std::numbers::pi_v<W>/2&&ps<std::numbers::pi_v<W>,"independent signed closed phase");
  compare(closed.rows[0],dh*std::sin(pl)/q/1.5L,dh*std::sin(ps)/q/101,
          dh*std::sin(ps-pl)/q/101);
  need(*closed.rows[0].outputs[3].value>1,"curved beta greater than one is physical");
  const LensEpochPair antipode{.5,10000};
  need(ThinLensGeometry(CurvedFLRW({70,.3,0,1.4})).evaluate({&antipode,1},all).status==S::outside_domain,
       "first conjugate point refused");
  const LensEpochPair moderate{.5,2};
  const auto flat=ThinLensGeometry(CurvedFLRW({70,.3,0,.7})).evaluate({&moderate,1},all);
  need(flat.status==S::ok,"flat control");
  for(double epsilon:{-1e-12,1e-12}) {
    const auto near=ThinLensGeometry(CurvedFLRW({70,.3,0,.7+epsilon})).evaluate({&moderate,1},all);
    need(near.status==S::ok,"near-flat continuity admission");
    for(unsigned n=0;n<6;++n)need(std::abs(*near.rows[0].outputs[n].value-*flat.rows[0].outputs[n].value)<=
        1e-9*std::abs(*flat.rows[0].outputs[n].value),"near-flat series both signs");
  }
  auto state=prepare_effective_neutrino_state({67.4,.02237,.12,2.7255,3.046,{{.06,.71611,1}}});
  auto request=state.observable_request(1059.95,"synthetic supplied drag","declared effective FD");
  need(request.has_value(),"thermal source request");
  auto provider=prepare_thermal_observables(*request);
  ThinLensGeometry thermal(provider);
  need(thermal.status()==S::ok&&thermal.thermal_provider()&&!thermal.curved_provider()&&
       thermal.provider_id()!=ThinLensGeometry(CurvedFLRW({70,.3,0,.7})).provider_id(),"distinct source identities");
  provider=ThermalObservables{};request.reset();
  const auto t=thermal.evaluate({&moderate,1},all);
  need(t.status==S::ok&&t.rows[0].momentum_callbacks>0,"retained FD consumer after source destruction");
  need(t.callbacks==t.outer_callbacks+t.momentum_callbacks&&t.callbacks==t.rows[0].callbacks,"FD counted once");
  const auto only=thermal.evaluate({&moderate,1},lens_geometry_mask(LensGeometryOutput::critical_surface_density_kg_m2));
  need(only.status==S::ok&&only.rows[0].outputs[4].value,"requested Sigma only");
  for(unsigned i:{0u,1u,2u,3u,5u})need(!only.rows[0].outputs[i].value&&only.rows[0].outputs[i].status==S::invalid_input,"omitted outputs absent");
  ThinLensGeometry copy=thermal,moved=std::move(copy);
  need(copy.status()==S::invalid_input&&moved.status()==S::ok,"move source invalidation");
  auto &alias=moved;moved=std::move(alias);need(moved.status()==S::ok,"self move");
  ThinLensGeometry simple(CurvedFLRW({70,.3,0,.7}));
  const LensEpochPair invalid[]{{0,1},{1,1},{2,1},{-1,2},{.5,std::numeric_limits<double>::quiet_NaN()},{.5,10001}};
  auto bad=simple.evaluate(invalid,all);need(bad.rows.size()==6&&bad.status!=S::ok,"all invalid rows retained");
  for(auto &r:bad.rows)for(auto &v:r.outputs)need(!v.value&&v.status!=S::ok,"refused outputs withheld");
  need(simple.evaluate({&moderate,1},0).rows.empty()&&simple.evaluate({&moderate,1},all+1).rows.empty(),"invalid masks");
  LensGeometryPolicy p;p.maximum_pairs=0;
  need(simple.evaluate({&moderate,1},all,p).status==S::work_limit,"pair quota");
  p={};p.maximum_native_bytes=1;need(simple.evaluate({&moderate,1},all,p).status==S::work_limit,"payload quota");
  need(ThinLensGeometry(*thermal.thermal_provider(),p).status()==S::work_limit,"retained provider payload quota");
  p={};p.maximum_callbacks_per_pair=3;
  const auto failed=simple.evaluate({&moderate,1},all,p);
  need(failed.status==S::work_limit&&failed.callbacks>0&&failed.callbacks<=3,"failed callbacks retained original quota");
  p={};p.maximum_total_callbacks=3;
  const auto total=simple.evaluate(pairs,all,p);need(total.rows.size()==4&&total.callbacks<=3&&total.status==S::work_limit,"whole batch quota");
  p={};p.relative_tolerance=1e-18;p.absolute_tolerance_mpc=1e-25;
  need(simple.evaluate({&moderate,1},all,p).status!=S::ok,"unresolved arithmetic refused");
  const auto before=std::fegetround();std::fesetround(FE_DOWNWARD);
  need(simple.evaluate({&moderate,1},all).status!=S::ok,"unsupported rounding");std::fesetround(before);
  std::cout<<"PASS thin lens geometry "<<checks<<" checks\n";
}
