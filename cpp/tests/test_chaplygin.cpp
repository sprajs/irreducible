// Continuity-derived fluid density ODE, independent Gauss-Legendre distance,
// exact endpoints and direct-algebra 100/160-digit binary64-input fixtures.
#include "irred/chaplygin.hpp"
#include "irred/quantities.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
using namespace irred::cosmology;
using W=long double;
using S=irred::numerics::Status;
namespace {
unsigned checks=0;
void need(bool b,const char *s){++checks;if(!b)throw std::runtime_error(s);}
void near(W a,W b,W tol,const char *s){need(std::abs(a-b)<=tol*std::max(1.0L,std::abs(b)),s);}
void relative(W a,W b,W tol,const char *s){need(std::abs(a-b)<=tol*std::abs(b),s);}
W direct_f(const ChaplyginModel &m,W a){return std::pow(W(m.a_s)+(1-W(m.a_s))*std::pow(a,-3*(1+W(m.alpha))),1/(1+W(m.alpha)));}
W direct_e(const ChaplyginModel &m,W a){
 const W cg=1-(W(m.ordinary_matter_fraction)+m.radiation_fraction);
 return std::sqrt(m.ordinary_matter_fraction/std::pow(a,3)+m.radiation_fraction/std::pow(a,4)+cg*direct_f(m,a));
}
// F_x=-3(F-A_s F^-alpha), integrated from F(1)=1. Production is closed algebra.
W density_ode(const ChaplyginModel &m,W a,unsigned steps){
 const W h=std::log(a)/steps;W f=1;
 auto rhs=[&](W y){return -3*(y-m.a_s*std::pow(y,-W(m.alpha)));};
 for(unsigned i=0;i<steps;++i){auto k1=rhs(f),k2=rhs(f+h*k1/2),k3=rhs(f+h*k2/2),k4=rhs(f+h*k3);f+=h*(k1+2*k2+2*k3+k4)/6;}
 return f;
}
// Gauss4 nodes/weights derived from P4=(35x^4-30x^2+3)/8, composite log-scale.
W gauss_distance(const ChaplyginModel &m,W a,unsigned cells){
 const W x0=std::log(a),h=-x0/cells,q=std::sqrt(6.0L/5);
 const W lo=std::sqrt((3-2*q)/7),hi=std::sqrt((3+2*q)/7);
 const W wl=(18+std::sqrt(30.0L))/36,wh=(18-std::sqrt(30.0L))/36;
 const std::array<W,4> nodes{-hi,-lo,lo,hi},weights{wh,wl,wl,wh};W sum=0;
 for(unsigned k=0;k<cells;++k)for(unsigned j=0;j<4;++j){
  const W x=x0+(k+.5L)*h+nodes[j]*h/2,b=std::exp(x);
  sum+=h/2*weights[j]/(b*direct_e(m,b));
 }
 return sum;
}
}
int main(){
 try{
  ChaplyginModel m;const W unit=(W(irred::speed_of_light_m_per_s)/1000)/m.anchor_hubble_km_s_mpc;
  const double scales[]{1,.8,.3,1e-4};
  for(double alpha:{0.,.4,1.})for(double as:{0.,1.}){
   m.a_s=as;m.alpha=alpha;
   auto out=evolve_chaplygin(m,scales,chaplygin_distances);need(out.status==S::ok,"analytic endpoints accepted");
   for(const auto &r:out.rows){
    const W a=r.scale_factor;
    relative(r.expansion_over_anchor_hubble,as==0?std::pow(a,-1.5L):1,3e-15L,"endpoint expansion");
    need(r.fluid_equation_of_state==-as && r.fluid_barotropic_slope==alpha*as,"endpoint EOS/slope");
    near(r.deceleration,as==0?.5L:-1,3e-15L,"endpoint deceleration");
    need(r.distances.has_value(),"requested endpoint distances");
    const W chi=as==0?2*(1-std::sqrt(a)):1/a-1;
    near(r.distances->comoving_mpc/unit,chi,3e-9L,"analytic endpoint comoving distance");
    near(r.distances->luminosity_mpc*a,r.distances->comoving_mpc,3e-15L,"flat luminosity relation");
    near(r.distances->angular_diameter_mpc/a,r.distances->comoving_mpc,3e-15L,"flat angular relation");
   }
  }
  // alpha=0 is exactly a density sum of vacuum and pressureless CG components,
  // in addition to the separately supplied ordinary matter and radiation.
  m={70,.2,.05,.73,0};auto alpha_zero=evolve_chaplygin(m,scales);
  need(alpha_zero.status==S::ok && alpha_zero.callbacks==0,"alpha0 background only");
  for(const auto &r:alpha_zero.rows){
   const W a=r.scale_factor,cg=alpha_zero.anchor_fluid_fraction;
   const W reference=(W(m.ordinary_matter_fraction)+cg*(1-W(m.a_s)))/std::pow(a,3)+m.radiation_fraction/std::pow(a,4)+cg*m.a_s;
   relative(r.expansion_over_anchor_hubble*r.expansion_over_anchor_hubble,reference,3e-15L,"alpha0 Lambda+dust equivalence");
   need(r.fluid_barotropic_slope==0 && !r.distances,"zero slope and unrequested distance absent");
  }
  const double mixed[]{1,.7,.7,.2,.02};
  for(double alpha:{0.,.5,1.})for(double as:{.01,.7,.999}){
   m={70,.2,.05,as,alpha};auto out=evolve_chaplygin(m,mixed,chaplygin_distances);
   need(out.status==S::ok,"mixed Chaplygin accepted");
   ChaplyginPolicy p;p.absolute_distance_tolerance_mpc=1e-9;p.relative_tolerance=1e-11;
   auto refined=evolve_chaplygin(m,mixed,chaplygin_distances,p);need(refined.status==S::ok,"native distance refinement");
   for(std::size_t k=0;k<out.rows.size();++k){
    const auto &r=out.rows[k];const W a=r.scale_factor;
    const W density=density_ode(m,a,8192),density2=density_ode(m,a,16384);
    // Frozen 2e-8 comparison, each independent refinement <=5% (1e-9).
    relative(density,density2,1e-9L,"independent density ODE refinement");
    relative(r.fluid_density_over_anchor_fluid_density,density2,2e-8L,"independent density ODE");
    const W chi=gauss_distance(m,a,128),chi2=gauss_distance(m,a,256);
    near(chi,chi2,1e-9L,"independent Gauss distance refinement");
    near(r.distances->comoving_mpc/unit,chi2,2e-8L,"independent Gauss distance");
    near(r.distances->comoving_mpc/unit,refined.rows[k].distances->comoving_mpc/unit,3e-9L,"native distance refinement comparison");
    near(r.ordinary_matter_fraction+r.radiation_fraction+r.fluid_fraction,1,3e-15L,"composition closure");
    relative(-r.fluid_equation_of_state,W(as)/std::pow(r.fluid_density_over_anchor_fluid_density,1+W(alpha)),3e-14L,"barotropic EOS identity");
    near(r.fluid_one_plus_equation_of_state-r.fluid_equation_of_state,1,3e-15L,"stable enthalpy decomposition");
    near(r.fluid_barotropic_slope,-alpha*r.fluid_equation_of_state,3e-15L,"formal EOS derivative");
    need(r.fluid_fraction>0 && r.fluid_equation_of_state<0 && r.fluid_one_plus_equation_of_state>0,"positive fluid and enthalpy");
   }
   need(out.rows[1].distances->comoving_mpc==out.rows[2].distances->comoving_mpc,"duplicate epoch retained distance");
  }
  // mpmath1.3.0, 100/160 dps, exact binary64 inputs via integer ratios,
  // direct power of the two-term sum (not production log-sum-exp). Reference
  // refinement was below1e-90 scaled. Literals retain >binary64 significance.
  struct Fixture{double as,alpha,a;W f,w,one,e;};
  const Fixture fixtures[]{
   {0x1p-52,1,1e-4,999999999999.999745212489465906481585374L,-2.22044604925031421233110566368303202384e-40L,0.999999999999999999999999999999999999999778L,22381912.3401017709304788113228663592449L},
   {0x1.fffffffffffffp-1,1,.001,10.5840588841198176728053296270641175546L,-.00892679384388314141722378949331427668086L,.991073206156116858582776210506685723319L,224053.565041795402300732076929937485861L},
   {0x1.fffffffffffffp-1,1,.01,1.00005550961057276933868382819620932453L,-.999888990022120824806458582029485549633L,.000111009977879175193541417970514450366573L,2280.35101465577708277631477928118751539L},
   {.7,.37,.003,15380587.2663943996317464920457928982858L,-9.97574047427976399115763899488723718517e-11L,.999999999900242595257202360088423610051L,25223.5365972832433883042664277709482786L}
  };
  for(const auto &f:fixtures){
   m={70,.2,.05,f.as,f.alpha};const double a[]{f.a};auto out=evolve_chaplygin(m,a);need(out.status==S::ok,"critical As fixture accepted");
   const auto &r=out.rows.front();
   relative(r.fluid_density_over_anchor_fluid_density,f.f,3e-14L,"high precision critical F");
   relative(r.fluid_equation_of_state,f.w,3e-14L,"high precision critical w");
   relative(r.fluid_one_plus_equation_of_state,f.one,3e-14L,"high precision critical enthalpy");
   relative(r.expansion_over_anchor_hubble,f.e,3e-14L,"high precision critical E");
  }
  m={70,.2,.05,.7,.5};auto kept=evolve_chaplygin(m,mixed,chaplygin_distances);auto copy=kept;kept.rows.clear();need(copy.rows.size()==5,"retained copy lifetime");
  auto moved=std::move(copy);need(moved.rows.size()==5,"retained move lifetime");
  need(evolve_chaplygin(m,{}).status==S::ok,"empty background batch");
  ChaplyginPolicy cap;cap.maximum_callbacks=0;need(evolve_chaplygin(m,mixed,0,cap).status==S::ok,"background only requires no callbacks");
  cap.maximum_callbacks=4;auto failed=evolve_chaplygin(m,mixed,chaplygin_distances,cap);
  need(failed.status==S::work_limit && failed.rows.empty() && failed.callbacks==4,"callback cap atomic refusal");
  cap={};cap.maximum_depth=0;failed=evolve_chaplygin(m,mixed,chaplygin_distances,cap);need(failed.status==S::work_limit && failed.rows.empty(),"depth cap atomic refusal");
  cap={};cap.maximum_native_bytes=1;need(evolve_chaplygin(m,mixed,0,cap).status==S::work_limit,"payload cap");
  cap={};cap.maximum_rows=2;need(evolve_chaplygin(m,mixed,0,cap).status==S::work_limit,"row cap");
  cap={};cap.relative_tolerance=0;need(evolve_chaplygin(m,mixed,0,cap).status==S::invalid_input,"precision admission");
  need(evolve_chaplygin(m,mixed,2).status==S::invalid_input,"output mask admission");
  auto bad=m;bad.ordinary_matter_fraction=.95;bad.radiation_fraction=.05;
  const auto tiny=evolve_chaplygin(bad,mixed);need(tiny.status==S::ok && tiny.anchor_fluid_fraction>0,"represented positive CG is never deleted");
  bad.ordinary_matter_fraction=.75;bad.radiation_fraction=.25;need(evolve_chaplygin(bad,mixed).status==S::outside_domain,"positive CG component required");
  bad=m;bad.radiation_fraction=-.01;need(evolve_chaplygin(bad,mixed).status==S::outside_domain,"negative ordinary radiation");
  bad=m;bad.a_s=1.01;need(evolve_chaplygin(bad,mixed).status==S::outside_domain,"As domain");
  bad=m;bad.alpha=-.1;need(evolve_chaplygin(bad,mixed).status==S::outside_domain,"alpha domain");
  bad=m;bad.anchor_hubble_km_s_mpc=0;need(evolve_chaplygin(bad,mixed).status==S::outside_domain,"physical H units/domain");
  bad=m;bad.a_s=std::numeric_limits<double>::quiet_NaN();need(evolve_chaplygin(bad,mixed).status==S::nonfinite_input,"nonfinite model");
  const double past[]{1,.5,.6};need(evolve_chaplygin(m,past).status==S::outside_domain,"past order");
  const double early[]{1e-5};need(evolve_chaplygin(m,early).status==S::outside_domain,"early cap");
  const double future[]{1.1};need(evolve_chaplygin(m,future).status==S::outside_domain,"future cap");
  const double nan[]{std::numeric_limits<double>::infinity()};need(evolve_chaplygin(m,nan).status==S::nonfinite_input,"nonfinite query");
  const int rounding=std::fegetround();std::fesetround(FE_UPWARD);const auto round_fail=evolve_chaplygin(m,mixed);std::fesetround(rounding);need(round_fail.status==S::invalid_input,"rounding profile");
  std::cout<<checks<<" Chaplygin scientific/resource checks passed\n";
 }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
