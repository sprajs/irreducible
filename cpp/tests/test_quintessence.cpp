#include "irred/quintessence.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using S=irred::numerics::Status;
namespace {
unsigned checks=0;
void need(bool b,const char *s){++checks;if(!b)throw std::runtime_error(s);}
void near(double a,long double b,double tol,const char *s){need(std::abs(a-b)<=tol*(1+std::abs(b)),s);}
// Independent dimensional Klein-Gordon/Friedmann variables, not autonomous
// density fractions: q=phi_dot/(Mpl Hanchor), U=V/(3 Mpl² Hanchor²).
std::array<long double,3> field_reference(ExponentialQuintessence p,double a,unsigned count) {
 using W=long double;W x=p.signed_kinetic_fraction_root,m=1-x*x-p.potential_fraction-p.radiation_fraction;
 std::array<W,3> y{std::sqrt(6.L)*x,p.potential_fraction,0};
 W h=std::log(W(a))/count;
 auto rhs=[&](W n,std::array<W,3> z){W e=std::sqrt(m*std::exp(-3*n)+p.radiation_fraction*std::exp(-4*n)+z[0]*z[0]/6+z[1]);return std::array<W,3>{-3*z[0]+3*p.lambda*z[1]/e,-p.lambda*z[1]*z[0]/e,std::exp(-n)/e};};
 auto add=[](std::array<W,3> y,std::array<W,3> k,W h){for(unsigned i=0;i<3;++i)y[i]+=h*k[i];return y;};
 for(unsigned j=0;j<count;++j){W n=j*h;auto k1=rhs(n,y),k2=rhs(n+h/2,add(y,k1,h/2)),k3=rhs(n+h/2,add(y,k2,h/2)),k4=rhs(n+h,add(y,k3,h));for(unsigned i=0;i<3;++i)y[i]+=h*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;}
 W n=std::log(W(a)),e=std::sqrt(m*std::exp(-3*n)+p.radiation_fraction*std::exp(-4*n)+y[0]*y[0]/6+y[1]);
 return {e,(y[0]*y[0]/6-y[1])/(y[0]*y[0]/6+y[1]),-299792.458L/p.h_anchor_km_s_mpc*y[2]};
}
}
int main(){try {
 const double axes[]{1,.7,.2,1.5};
 // Frozen comparison allocations: analytic3e-8, independent transient2e-7;
 // independent 4096→8192 refinement <=1% of transient allocation.
 for(auto p:{ExponentialQuintessence{0,70,0,.7,.0001},ExponentialQuintessence{0,70,1,0,0},ExponentialQuintessence{0,70,-1,0,0},ExponentialQuintessence{std::sqrt(1.5),70,.5,.75,0},ExponentialQuintessence{std::sqrt(6.),70,.5,.25,0}}) {
  auto owner=prepare_quintessence(p);need(owner.status()==S::ok,"analytic initial admission");auto b=owner.evaluate(axes,true);need(b.status==S::ok&&b.rows.size()==4,"analytic ordered batch");
  for(size_t i=0;i<4;++i){auto &r=b.rows[i];need(r.e.value&&r.h_km_s_mpc.value&&r.omega_phi.value,"background available");long double a=axes[i],e=0,w=0;
   if(p.lambda==0&&p.potential_fraction>.5){e=std::sqrt(.2999L/std::pow(a,3)+.0001L/std::pow(a,4)+.7L);w=-1;}
   else if(std::abs(p.signed_kinetic_fraction_root)==1){e=1/std::pow(a,3);w=1;}
   else if(p.potential_fraction==.75){e=std::pow(a,-.75L);w=-.5;}
   else {e=std::pow(a,-1.5L);w=0;}
   near(*r.e.value,e,3e-8,"analytic E");near(*r.h_km_s_mpc.value,70*e,3e-8,"analytic H");near(*r.w_phi.value,w,3e-8,"analytic w");need(r.constraint_residual<=1.1e-9,"Friedmann constraint");
   if(std::abs(p.signed_kinetic_fraction_root)==1&&a<=1)near(*r.dm_mpc.value,299792.458L/70*(1-a*a)/2,3e-8,"stiff exact distance");
   if(a>1)need(r.dm_mpc.status==S::outside_domain&&r.e.status==S::ok,"future distance retains background");
  }
 }
 // Radiation fixed point at lambda4; represented input closure retained.
 ExponentialQuintessence radiation{4,70,std::sqrt(1./6),1./12,.75};
 auto rb=prepare_quintessence(radiation).evaluate(std::span(axes).first(3));need(rb.status==S::ok,"radiation scaling batch");for(size_t i=0;i<3;++i){near(*rb.rows[i].e.value,1.L/(axes[i]*axes[i]),3e-8,"radiation scaling E");near(*rb.rows[i].w_phi.value,1.L/3,3e-8,"radiation scaling w");}
 ExponentialQuintessence transient{1.2,67.4,.08,.65,.00009};auto owner=prepare_quintessence(transient);auto b=owner.evaluate(axes,true);need(b.status==S::ok,"transient batch");
 for(size_t i=0;i<4;++i){auto ref=field_reference(transient,axes[i],8192),coarse=field_reference(transient,axes[i],4096);for(size_t k=0;k<3;++k)need(std::abs(ref[k]-coarse[k])<2e-9*(1+std::abs(ref[k])),"independent refinement");near(*b.rows[i].e.value,ref[0],2e-7,"independent E");near(*b.rows[i].w_phi.value,ref[1],2e-7,"independent w");if(axes[i]<=1)near(*b.rows[i].dm_mpc.value,ref[2],2e-7,"independent distance");}
 // Independent composite Simpson quadrature for the lambda0 distance.
 auto lcdm=prepare_quintessence({0,70,0,.7,.0001}).evaluate(std::span(axes).subspan(2,1),true);
 long double integral=0,low=axes[2],step=(1-low)/16384;
 for(unsigned j=0;j<=16384;++j){long double a=low+j*step,e=std::sqrt(.2999L/(a*a*a)+.0001L/(a*a*a*a)+.7L);integral+=(j==0||j==16384?1:j%2?4:2)/(a*a*e);}
 near(*lcdm.rows[0].dm_mpc.value,299792.458L/70*integral*step/3,3e-8,"independent Lambda distance quadrature");
 auto reflected=transient;reflected.lambda=-transient.lambda;reflected.signed_kinetic_fraction_root=-transient.signed_kinetic_fraction_root;
 auto reflection=prepare_quintessence(reflected).evaluate(axes,true);need(reflection.status==S::ok,"field reflection admission");for(size_t i=0;i<4;++i)near(*reflection.rows[i].e.value,*b.rows[i].e.value,1e-12,"lambda phi reflection");
 auto tight=QuintessencePolicy{};tight.absolute_tolerance=1e-12;tight.relative_tolerance=1e-11;auto refined=owner.evaluate(axes,true,tight);need(refined.status==S::ok,"producer refinement");for(size_t i=0;i<4;++i)near(*b.rows[i].e.value,*refined.rows[i].e.value,2e-7,"producer E refinement");
 auto zero=prepare_quintessence({0,70,0,0,0}).evaluate(std::span(axes).first(1));need(zero.rows[0].w_phi.status==S::singular&&zero.rows[0].e.value,"zero scalar preserves E");
 for(auto invalid:{ExponentialQuintessence{21},ExponentialQuintessence{0,0},ExponentialQuintessence{0,70,1,.1},ExponentialQuintessence{0,70,0,-.1}})need(prepare_quintessence(invalid).status()==S::outside_domain,"invalid physical state");
 auto nan=transient;nan.lambda=std::numeric_limits<double>::quiet_NaN();need(prepare_quintessence(nan).status()==S::nonfinite_input,"nonfinite state");
 auto limited=QuintessencePolicy{};limited.maximum_total_callbacks=12;auto fail=owner.evaluate(axes,true,limited);need(fail.status==S::work_limit&&fail.callbacks<=12&&fail.rows.size()==4,"all rows retained at cap");need(fail.rows[0].e.value&&fail.rows[1].e.status==S::work_limit,"partial rows preserved");
 limited=QuintessencePolicy{};limited.maximum_native_bytes=1;need(owner.evaluate(axes,false,limited).status==S::work_limit,"payload refusal");need(!quintessence_payload_bound(std::numeric_limits<size_t>::max()),"payload overflow");
 const double bad[]{0,std::numeric_limits<double>::quiet_NaN(),1};auto mixed=owner.evaluate(bad);need(mixed.rows[0].status==S::outside_domain&&mixed.rows[1].status==S::nonfinite_input&&mixed.rows[2].e.value,"invalid axes preserve valid row");
 auto impossible=QuintessencePolicy{};impossible.absolute_tolerance=1e-30;impossible.relative_tolerance=0;impossible.maximum_halvings=2;need(owner.evaluate(std::span(axes).subspan(1,1),false,impossible).status==S::conditioning_budget_exceeded,"arithmetic impossible policy refused");
 std::fesetround(FE_DOWNWARD);need(owner.evaluate(axes).status==S::invalid_input,"rounding refusal");std::fesetround(FE_TONEAREST);
 std::cout<<checks<<" quintessence checks passed\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
