#include "hydrogen_thermal_reference.hpp"
#ifndef IRRED_REFERENCE_ONLY
#include "irred/recombination_drag.hpp"
#endif
#include <iomanip>
#include <iostream>
#include <string>
namespace {
namespace ref=hydrogen_thermal_reference;
using ref::W;
unsigned checks=0;
void require(bool value,const std::string& what) {
 ++checks;if(!value)throw std::runtime_error(what);
}
W budget(unsigned group,W value) {
 constexpr W absolute[]{1e-8L,2e-7L,2e-5L,2e-7L,2e-9L,2e-9L,0};
 constexpr W relative[]{1e-6L,1e-6L,2e-7L,1e-6L,2e-6L,3e-6L,1e-5L};
 return absolute[group]+relative[group]*std::abs(value);
}
std::array<W,7> coordinates(ref::Value a) {
 return {a.x,a.drag,a.temperature,a.thomson,a.opacity,a.visibility,a.survival};
}
void analytic() {
 const W nonfinite=std::numeric_limits<W>::quiet_NaN();
 require(!std::isfinite(ref::norm({0,0,nonfinite,0})),"nonfinite residual norm refuses rigorous zero");
 bool refused=false;ref::Stats invalid_stats;
 try {
  ref::radau([=](W,const ref::State&){return ref::State{nonfinite,0};},0,.01L,{.5L,1},invalid_stats);
 } catch(const std::runtime_error&) {refused=true;}
 require(refused&&invalid_stats.steps==0,"nonfinite stage equation explicitly refused");
 for(W A:{0.L,.5L,10.L,1e6L}) {
  ref::Stats stats;ref::State y{.5L,1};W span=A==0 ? 1 : 12/(1+A);
  unsigned steps=A==0 ? 8192 : 2048;W ds=span/steps;
  auto equation=[A](W,const ref::State& y){return ref::State{0,A-(1+A)*y[1]};};
  for(unsigned i=0;i<steps;++i)y=ref::radau(equation,i*ds,ds,y,stats).end;
  W exact=A/(1+A)+(1-A/(1+A))*std::exp(-(1+A)*span);
  require(std::abs(y[1]/exact-1)<2e-12L,"exact constant Compton temperature");
  require(stats.min_x>0&&stats.min_neutral>0&&stats.min_theta>0,"analytic stage positivity");
  require(stats.residual<=8e-17L,"analytic converged stage residual");
  if(A==0)require(std::abs(y[1]*std::exp(-span)/std::exp(-2*span)-1)<2e-12L,"adiabatic Kelvin scaling");
 }
 // Original redshift-density controls, outside the physical source admission.
 // These check the quadrature/finite-mass algebra, not a new cosmology.
 for(int power:{-1,0,2}) {
  W lo=2,hi=5,k=.02L;
  auto depth=[=](W z){return power==-1 ? k*std::log((1+z)/(1+lo)) : k/(power+1)*(std::pow(1+z,power+1)-std::pow(1+lo,power+1));};
  W total=depth(hi),mass=0;
  for(unsigned i=0;i<16;++i)mass+=ref::gl8([&](W z){return k*std::pow(1+z,power)*std::exp(-depth(z));},lo+(hi-lo)*i/16,lo+(hi-lo)*(i+1)/16);
  require(std::abs(mass/(-std::expm1(-total))-1)<2e-12L,"analytic opacity visibility integral");
  require(std::abs(mass+std::exp(-total)-1)<2e-12L,"visible plus surviving finite mass");
  require(depth(lo)==0,"analytic supplied endpoint depth zero");
 }
}
struct Case { const char* name;ref::Model model;double initial,late; };
void physical(const Case& test) {
 std::vector<W> z{test.initial,test.initial-.001,test.initial-.01,test.initial-.1};
 for(double value:{1500.,1400.,1300.,1200.,1100.,1000.,900.,800.,700.,600.,500.,400.})
  if(value<test.initial-.1 && value>test.late)z.push_back(value);
 z.push_back(test.late);
 auto coarse=ref::integrate(test.model,test.initial,test.late,z,1);
 auto medium=ref::integrate(test.model,test.initial,test.late,z,2);
 auto fine=ref::integrate(test.model,test.initial,test.late,z,4);
 W maximum_fraction=0,initial_fraction=0;
 for(size_t i=0;i<z.size();++i) {
  auto a=coordinates(medium.query(i)),b=coordinates(fine.query(i));
  for(unsigned group=0;group<7;++group) {
   W delta=std::abs(a[group]-b[group]),allowance=budget(group,b[group]);
   W fraction=allowance>0 ? delta/allowance : delta==0 ? 0 : std::numeric_limits<W>::infinity();
   maximum_fraction=std::max(maximum_fraction,fraction);if(i<4)initial_fraction=std::max(initial_fraction,fraction);
   if(!(fraction<=.05L))std::cerr<<"REFINEMENT case="<<test.name<<" z="<<z[i]<<" group="<<group<<" medium="<<a[group]<<" fine="<<b[group]<<" allocation_fraction="<<fraction<<'\n';
   require(fraction<=.05L,"reference refinement within5% downstream allocation");
  }
  auto q=fine.query(i);
  require(q.x>0&&q.x<1&&q.temperature>0&&q.survival>0&&q.visibility>0,"reference output positivity");
  require(q.temperature<=test.model.T0*(1+z[i])*(1+1e-15L),"reference physical temperature below radiation");
 }
 for(const auto* r:{&coarse,&medium,&fine}) {
  require(r->stats.residual<=8e-17L,"physical converged stage residual");
  require(r->stats.min_x>0&&r->stats.min_neutral>0&&r->stats.min_theta>0,"physical stage positivity");
  require(r->stats.maximum_initial_step_stiffness<=.125L*(1+1e-15L),"initial Compton layer explicitly resolved");
 }
 require(std::abs(medium.root-fine.root)<=.05L*.002L,"independent drag-root refinement");
 auto quad=ref::opacity_gl(fine);
 require(std::abs(quad.thomson-fine.nodes.back().thomson)<=.05L*budget(3,quad.thomson),"independent GL Thomson depth");
 require(std::abs(quad.drag-fine.nodes.back().drag)<=.05L*budget(1,quad.drag),"independent GL drag depth");
 require(std::abs(quad.survival-fine.query(0).survival)<=.05L*budget(6,quad.survival),"independent GL positive survivor");
 require(std::abs(quad.mass+quad.survival-1)<=.05L*2e-7L,"finite-interval GL mass plus survivor");
 std::cout<<"REFERENCE "<<test.name<<" steps="<<coarse.stats.steps<<','<<medium.stats.steps<<','<<fine.stats.steps
  <<" rhs="<<fine.stats.rhs<<" Newton="<<fine.stats.newton<<" residual="<<fine.stats.residual
  <<" minimum_stage_x="<<fine.stats.min_x<<" minimum_stage_neutral="<<fine.stats.min_neutral
  <<" minimum_stage_theta="<<fine.stats.min_theta<<" initial_A_bound="<<fine.stats.initial_stiffness
  <<" initial_hA="<<fine.stats.maximum_initial_step_stiffness<<" max_refinement_fraction="<<maximum_fraction
  <<" boundary_refinement_fraction="<<initial_fraction<<" root_delta="<<std::abs(medium.root-fine.root)
  <<" GL_tau_delta="<<std::abs(quad.thomson-fine.nodes.back().thomson)
  <<" GL_drag_delta="<<std::abs(quad.drag-fine.nodes.back().drag)
  <<" GL_mass_error="<<std::abs(quad.mass+quad.survival-1)<<" survivor="<<quad.survival<<'\n';
#ifndef IRRED_REFERENCE_ONLY
 using namespace irred::cosmology;
 ThermalPhysicalModel model{double(test.model.H0),double(test.model.baryon),double(test.model.cdm),double(test.model.T0),double(test.model.other),{}};
 for(auto a:test.model.species)model.species.push_back({0,double(a.temperature),double(a.weight)});
 PureHydrogenRequest request{model,test.initial,test.late,HydrogenTemperatureModel::evolved_compton_adiabatic};
 auto owner=prepare_pure_hydrogen_history(request);
 require(owner.status()==irred::numerics::Status::ok,"evolved native preparation");
 std::vector<double> queries;for(W value:z)queries.push_back(double(value));
 auto native=owner.evaluate(queries,127);
 require(native.status==irred::numerics::Status::ok&&native.rows.size()==z.size(),"native ordered coarse batch");
 W maximum_native_fraction=0;
 for(size_t i=0;i<z.size();++i) {
  auto reference=coordinates(fine.query(i));const auto& row=native.rows[i];
  const HydrogenHistoryValue* groups[]{&row.electron_fraction,&row.drag_depth,&row.matter_temperature_kelvin,&row.thomson_depth,&row.thomson_opacity_per_redshift,&row.visibility_per_redshift,&row.finite_endpoint_survival};
  for(unsigned group=0;group<7;++group) {
   require(groups[group]->status==irred::numerics::Status::ok&&groups[group]->value.has_value(),"native requested group accepted");
   W fraction=std::abs(*groups[group]->value-reference[group])/budget(group,reference[group]);
   maximum_native_fraction=std::max(maximum_native_fraction,fraction);
   if(!(fraction<=1))std::cerr<<"NATIVE case="<<test.name<<" z="<<z[i]<<" group="<<group<<" value="<<*groups[group]->value<<" reference="<<reference[group]<<" allocation_fraction="<<fraction<<'\n';
   require(fraction<=1,"native independent Radau comparison");
  }
 }
 auto root=owner.conditional_unit_depth_redshift();
 require(root.status==irred::numerics::Status::ok&&root.value&&std::abs(*root.value-fine.root)<=.002L,"native conditional drag root");
 std::cout<<"NATIVE "<<test.name<<" maximum_allocation_fraction="<<maximum_native_fraction<<'\n';
 // Integrate actual native visibility, retaining the same prepared history.
 // Bounded batches request only g; the independent GL rule supplies weights.
 const auto prepared_work=owner.work().total();
 const W survivor=*native.rows.front().finite_endpoint_survival.value;
 ref::Gauss8 gauss;W previous_mass=0;
 for(unsigned panels:{128u,256u,512u}) {
  W mass=0,width=(test.initial-test.late)/W(panels);
  std::vector<double> points;std::vector<W> weights;
  points.reserve(1024);weights.reserve(1024);
  auto consume=[&] {
   auto values=owner.evaluate(points,hydrogen_visibility);
   require(values.status==irred::numerics::Status::ok&&values.rows.size()==points.size(),"bounded native visibility batch");
   for(size_t i=0;i<points.size();++i) {
    const auto& g=values.rows[i].visibility_per_redshift;
    require(g.status==irred::numerics::Status::ok&&g.value&&*g.value>0,"native quadrature visibility admitted positive");
    mass+=weights[i]* *g.value;
   }
   require(owner.work().total()==prepared_work,"visibility batches retain prepared evolution work");
   points.clear();weights.clear();
  };
  for(unsigned i=0;i<panels;++i)for(unsigned j=0;j<8;++j) {
   points.push_back(double(test.late+(i+.5L)*width+width/2*gauss.x[j]));
   weights.push_back(width/2*gauss.w[j]);
   if(points.size()==1024)consume();
  }
  if(!points.empty())consume();
  const W conservation=std::abs(mass+survivor-1);
  const W refinement=panels==128 ? 0 : std::abs(mass-previous_mass);
  std::cout<<"NATIVE_MASS "<<test.name<<" panels="<<panels<<" mass="<<mass<<" survivor="<<survivor
   <<" conservation_error="<<conservation<<" refinement="<<refinement<<'\n';
  require(conservation<=2e-7L,"actual native finite visible mass plus survivor");
  if(panels==512)require(refinement<=1e-8L,"actual native visibility GL refinement");
  previous_mass=mass;
 }
#endif
}
} // namespace
int main() {
 try {
  std::cout<<std::setprecision(21)<<std::unitbuf;std::cerr<<std::setprecision(21);analytic();
  ref::Model low{60.,.015,.08,2.70,0,{{1.9,2}}};
  ref::Model high{80.,.03,.15,2.75,3e-5,{}};
  physical({"fiducial",{},1600,300});physical({"low-density-massless-species",low,1550,600});
  physical({"high-density-late450",high,1600,450});physical({"high-density-late300",high,1600,300});
  std::cout<<"PASS "<<checks<<" independent hydrogen thermal controls; common atomic/SI assets and long-double/libm, original Radau/direct-SI/GL algorithms\n";
 } catch(const std::exception& e) {std::cerr<<"FAIL after "<<checks<<" controls: "<<e.what()<<'\n';return 1;}
}
