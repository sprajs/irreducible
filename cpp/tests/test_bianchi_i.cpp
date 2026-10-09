// Source-free Einstein/shear/null-momentum identities, exact Kasner controls,
// independent composite Gauss-Legendre quadrature and proper-time physical ODE.
#include "irred/bianchi_i.hpp"
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
void need(bool b,const char *s) { ++checks;if(!b)throw std::runtime_error(s); }
void near(W a,W b,W tol,const char *s) { need(std::abs(a-b)<=tol*std::max(1.0L,std::abs(b)),s); }
std::array<W,3> shear(const BianchiIModel &m) {
 return {m.shear_x_over_anchor_hubble,m.shear_y_over_anchor_hubble,
   -(W(m.shear_x_over_anchor_hubble)+m.shear_y_over_anchor_hubble)};
}
W sigma(const BianchiIModel &m) { auto s=shear(m);return (s[0]*s[0]+s[1]*s[1]+s[2]*s[2])/6; }
W expansion(const BianchiIModel &m,W a) {
 return std::sqrt(m.matter_fraction/(a*a*a)+m.radiation_fraction/std::pow(a,4)+
   m.lambda_fraction+sigma(m)/std::pow(a,6));
}
// Four-node Gauss-Legendre roots/weights derived from P4=(35x^4-30x^2+3)/8.
// Integrand changes with this reference only, not production's Simpson callback.
W gauss(const BianchiIModel &m,W a,unsigned segments,bool lookback=false) {
 const W x0=std::log(a),h=-x0/segments;
 const W q=std::sqrt(6.0L/5),lo=std::sqrt((3-2*q)/7),hi=std::sqrt((3+2*q)/7);
 const W wlo=(18+std::sqrt(30.0L))/36,whi=(18-std::sqrt(30.0L))/36;
 const std::array<W,4> nodes{-hi,-lo,lo,hi},weights{whi,wlo,wlo,whi};
 W sum=0;
 for(unsigned k=0;k<segments;++k) {
  const W mid=x0+(k+.5L)*h;
  for(unsigned j=0;j<4;++j) {
   const W b=std::exp(mid+nodes[j]*h/2),e=expansion(m,b);
   sum+=h/2*weights[j]*(lookback?1/e:1/(b*b*b*e));
  }
 }
 return sum;
}
using State=std::array<W,7>; // a, beta_i, dot(beta_i)/H_anchor
State proper_time(const BianchiIModel &m,W u,unsigned steps) {
 auto s=shear(m);State y{1,0,0,0,s[0],s[1],s[2]};
 const W h=u/steps;
 auto rhs=[&](State z) {
  const W a=z[0],sum=z[4]*z[4]+z[5]*z[5]+z[6]*z[6];
  const W e=std::sqrt(m.matter_fraction/std::pow(a,3)+m.radiation_fraction/std::pow(a,4)+m.lambda_fraction+sum/6);
  return State{-a*e,-z[4],-z[5],-z[6],3*e*z[4],3*e*z[5],3*e*z[6]};
 };
 for(unsigned i=0;i<steps;++i) {
  auto add=[&](State k,W f){State z{};for(unsigned j=0;j<7;++j)z[j]=y[j]+f*k[j];return z;};
  auto k1=rhs(y),k2=rhs(add(k1,h/2)),k3=rhs(add(k2,h/2)),k4=rhs(add(k3,h));
  for(unsigned j=0;j<7;++j)y[j]+=h*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;
 }
 return y;
}
}
int main() {
 try {
  BianchiIModel m; m.anchor_hubble_per_second=2e-18;
  const std::array<BianchiIQuery,4> q{{{1,{1,0,0}},{.8,{0,1,0}},{.4,{0,0,1}},{1e-4,{.6,.8,0}}}};
  // Exact isotropic FLRW limits, including anchor unit-frame direction.
  for(unsigned species=0;species<4;++species) {
   m.matter_fraction=species==0?1:species==3?.3:0;
   m.radiation_fraction=species==1?1:species==3?.1:0;
   m.lambda_fraction=species==2?1:species==3?.6:0;
   auto out=evolve_bianchi_i(m,q);need(out.status==S::ok,"isotropic model accepted");
   need(out.callbacks==0 && out.anchor_shear_fraction==0,"isotropic no unnecessary shear work");
   for(const auto &r:out.rows) {
    const W a=r.query.scale_factor;
    near(r.expansion_over_anchor_hubble,expansion(m,a),2e-15L,"isotropic FLRW expansion");
    near(r.one_plus_directional_redshift,1/a,2e-15L,"isotropic redshift");
    need(r.shear_fraction==0,"isotropic shear exactly absent");
    for(unsigned i=0;i<3;++i) {
     need(r.beta[i]==0,"isotropic beta exactly zero");
     near(r.directional_scale_factors[i],a,2e-15L,"isotropic axes");
     near(r.axis_hubble_per_second[i]/m.anchor_hubble_per_second,expansion(m,a),2e-15L,"isotropic axis H");
    }
   }
  }
  // Flat vacuum Kasner s=(2,-1,-1): physical axes (a^3,1,1).
  m.matter_fraction=m.radiation_fraction=m.lambda_fraction=0;
  m.shear_x_over_anchor_hubble=2;m.shear_y_over_anchor_hubble=-1;
  const std::array<BianchiIQuery,5> kasner{{{1,{1,0,0}},{.5,{1,0,0}},{.5,{0,1,0}},{.5,{.6,.8,0}},{1e-4,{1,0,0}}}};
  auto vac=evolve_bianchi_i(m,kasner);need(vac.status==S::ok && vac.callbacks==0,"Kasner analytic branch");
  for(const auto &r:vac.rows) {
   const W a=r.query.scale_factor;
   near(r.expansion_over_anchor_hubble,1/std::pow(a,3),3e-15L,"Kasner expansion");
   near(r.shear_fraction,1,3e-15L,"Kasner shear dominates");
   near(r.beta[0],2*std::log(a),3e-15L,"Kasner beta x");
   near(r.beta[1],-std::log(a),3e-15L,"Kasner beta y");
   near(r.directional_scale_factors[0],std::pow(a,3),3e-15L,"Kasner physical x axis");
   near(r.directional_scale_factors[1],1,3e-15L,"Kasner static transverse axis");
   near(r.axis_hubble_per_second[1]/r.mean_hubble_per_second,0,3e-15L,"Kasner zero axis H relative to mean expansion");
   const auto &n=r.query.observed_direction;
   near(r.one_plus_directional_redshift,std::sqrt(W(n[0])*n[0]/std::pow(a,6)+W(n[1])*n[1]+W(n[2])*n[2]),3e-15L,"Kasner axis/off-axis null momentum");
  }
  // Another vacuum Kasner s=(-2,1,1) has a contracting x axis: genuine blueshift.
  m.shear_x_over_anchor_hubble=-2;m.shear_y_over_anchor_hubble=1;
  auto contracting=evolve_bianchi_i(m,kasner);need(contracting.status==S::ok,"contracting axis accepted");
  need(contracting.rows[1].axis_hubble_per_second[0]<0 && contracting.rows[1].directional_redshift<0,"negative axis H and directional blueshift retained");
  for(const auto &r:contracting.rows) {
   const W a=r.query.scale_factor;
   const auto &n=r.query.observed_direction;
   near(r.one_plus_directional_redshift,std::sqrt(W(n[0])*n[0]*a*a+(W(n[1])*n[1]+W(n[2])*n[2])/std::pow(a,4)),3e-15L,"contracting Kasner null geometry");
  }
  // Named mixed-stress cases: production Simpson vs independently derived
  // Gauss4 and proper-time a/beta/shear ODE. Frozen 2e-8 comparison allocation;
  // each independent refinement must occupy <=5% (1e-9) of that allocation.
  const std::array<BianchiIQuery,5> mixed{{{1,{1,0,0}},{.7,{.6,.8,0}},{.7,{0,0,1}},{.4,{0,1,0}},{.3,{1,0,0}}}};
  for(W amplitude:{.01L,.2L,.8L}) {
   m.shear_x_over_anchor_hubble=static_cast<double>(amplitude);
   m.shear_y_over_anchor_hubble=static_cast<double>(-.3L*amplitude);
   m.matter_fraction=.25;m.radiation_fraction=.15;
   m.lambda_fraction=static_cast<double>(1-.4L-sigma(m));
   auto out=evolve_bianchi_i(m,mixed);need(out.status==S::ok,"mixed anisotropy accepted");
   BianchiIPolicy p;p.absolute_beta_tolerance=1e-12;p.relative_tolerance=1e-11;
   auto refined=evolve_bianchi_i(m,mixed,p);need(refined.status==S::ok,"native quadrature refinement");
   need(out.callbacks>0,"mixed callbacks accounted");
   for(std::size_t k=0;k<out.rows.size();++k) {
    const auto &r=out.rows[k];const W a=r.query.scale_factor;
    auto s=shear(m);
    const W integral=gauss(m,a,128),integral2=gauss(m,a,256);
    near(integral,integral2,1e-9L,"Gauss shear refinement");
    const W u=gauss(m,a,256,true),u2=gauss(m,a,512,true);
    near(u,u2,1e-9L,"Gauss elapsed time refinement");
    auto reference=proper_time(m,u2,8192),reference2=proper_time(m,u2,16384);
    for(unsigned j=0;j<7;++j)near(reference[j],reference2[j],1e-9L,"proper-time ODE refinement");
    near(reference2[0],a,2e-8L,"proper-time independent scale");
    near(r.shear_fraction,sigma(m)/std::pow(a,6)/std::pow(expansion(m,a),2),2e-15L,"physical shear ratio");
    near(r.directional_scale_factors[0]*r.directional_scale_factors[1]*r.directional_scale_factors[2],a*a*a,2e-13L,"trace-free physical volume");
    for(unsigned j=0;j<3;++j) {
     near(r.beta[j],-s[j]*integral2,2e-8L,"independent Gauss beta");
     near(r.beta[j],reference2[1+j],2e-8L,"proper-time independent beta");
     near(reference2[4+j],s[j]/std::pow(a,3),2e-8L,"proper-time shear conservation");
     near(r.beta[j],refined.rows[k].beta[j],3e-9L,"native beta refinement");
     near(r.axis_hubble_per_second[j]/m.anchor_hubble_per_second,expansion(m,a)+s[j]/std::pow(a,3),2e-15L,"axis expansion equation");
    }
    W energy2=0;
    for(unsigned j=0;j<3;++j)energy2+=W(r.query.observed_direction[j])*r.query.observed_direction[j]/std::pow(r.directional_scale_factors[j],2);
    near(r.one_plus_directional_redshift,std::sqrt(energy2),2e-15L,"geometric conserved photon momenta");
   }
   need(out.rows[1].beta==out.rows[2].beta,"duplicate epochs share anisotropy");
  }
  auto kept=evolve_bianchi_i(m,mixed);auto copy=kept;kept.rows.clear();need(copy.rows.size()==mixed.size(),"retained copy lifetime");
  auto moved=std::move(copy);need(moved.rows.size()==mixed.size(),"retained move lifetime");
  auto empty=evolve_bianchi_i(m,{});need(empty.status==S::ok && empty.callbacks==0,"empty batch");
  BianchiIPolicy cap;cap.maximum_callbacks=4;
  auto fail=evolve_bianchi_i(m,mixed,cap);need(fail.status==S::work_limit && fail.rows.empty() && fail.callbacks==4,"callback cap atomic refusal");
  cap={};cap.maximum_native_bytes=1;need(evolve_bianchi_i(m,mixed,cap).status==S::work_limit,"payload cap");
  cap={};cap.maximum_rows=2;need(evolve_bianchi_i(m,mixed,cap).status==S::work_limit,"row cap");
  cap={};cap.maximum_depth=0;
  auto depth_fail=evolve_bianchi_i(m,mixed,cap);need(depth_fail.status==S::work_limit && depth_fail.rows.empty(),"depth cap atomic refusal");
  cap={};cap.maximum_depth=61;need(evolve_bianchi_i(m,mixed,cap).status==S::invalid_input,"depth admission");
  cap={};cap.absolute_beta_tolerance=0;need(evolve_bianchi_i(m,mixed,cap).status==S::invalid_input,"precision admission");
  auto bad=m;bad.lambda_fraction+=.1;need(evolve_bianchi_i(bad,mixed).status==S::outside_domain,"no density renormalization");
  bad=m;bad.matter_fraction=-.1;need(evolve_bianchi_i(bad,mixed).status==S::outside_domain,"negative fluid refusal");
  bad=m;bad.anchor_hubble_per_second=0;need(evolve_bianchi_i(bad,mixed).status==S::outside_domain,"physical H required");
  bad=m;bad.shear_x_over_anchor_hubble=std::numeric_limits<double>::infinity();need(evolve_bianchi_i(bad,mixed).status==S::nonfinite_input,"nonfinite shear");
  auto invalid=mixed;invalid[0].observed_direction={1,1,0};need(evolve_bianchi_i(m,invalid).status==S::outside_domain,"nonunit direction no normalization");
  invalid=mixed;invalid[0].observed_direction={0,0,0};need(evolve_bianchi_i(m,invalid).status==S::outside_domain,"zero direction");
  invalid=mixed;invalid[0].observed_direction[0]=std::numeric_limits<double>::quiet_NaN();need(evolve_bianchi_i(m,invalid).status==S::nonfinite_input,"nonfinite direction");
  invalid=mixed;invalid[3].scale_factor=.8;need(evolve_bianchi_i(m,invalid).status==S::outside_domain,"past-ordered scale contract");
  invalid=mixed;invalid[0].scale_factor=1.1;need(evolve_bianchi_i(m,invalid).status==S::outside_domain,"observer endpoint domain");
  invalid=mixed;invalid[4].scale_factor=1e-5;need(evolve_bianchi_i(m,invalid).status==S::outside_domain,"early scale cap");
  bad=m;bad.shear_x_over_anchor_hubble=1e300;need(evolve_bianchi_i(bad,mixed).status==S::outside_domain,"huge shear closure refusal");
  const int rounding=std::fegetround();std::fesetround(FE_UPWARD);
  const auto rounding_result=evolve_bianchi_i(m,mixed);std::fesetround(rounding);
  need(rounding_result.status==S::invalid_input,"rounding profile admission");
  std::cout<<checks<<" Bianchi-I scientific/resource checks passed\n";
 }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
