// Independent CPL peer: derivatives of E^2, composite eight-point Gauss-Legendre
// in redshift (owner uses scale-factor Simpson), and constant-w exact limits.
// Shared conservation equations/libm are disclosed ancestry; not an interval
// proof. Fixed budgets from P01: integrals 2e-15+2e-10*abs(reference),
// analytic E/q/j 2e-13+2e-12*abs(reference), refinement <=1/10 budget.
#include "irred/background.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
unsigned checks=0; long double max_ratio=0,max_ref_ratio=0;
void check(bool x){++checks;if(!x)throw std::runtime_error("CPL peer check");}
void near(long double got,long double ref,long double abs,long double rel){
 auto budget=abs+rel*std::abs(ref);auto ratio=std::abs(got-ref)/budget;
 max_ratio=std::max(max_ratio,ratio);check(ratio<=1);
}
long double squared(long double z,long double om,long double w,long double wa){
 auto u=1+z;
 // Power and exponential separately, rather than production log exponential.
 return om*u*u*u+(1-om)*std::pow(u,3*(1+w+wa))*std::exp(-3*wa*z/u);
}
long double panels(long double z,long double om,long double w,long double wa,
                   unsigned count,bool clock){
 constexpr long double nodes[]{.183434642495649804939476142360L,.525532409916328985817739049189L,.796666477413626739591553936476L,.960289856497536231683560868569L};
 constexpr long double weights[]{.362683783378361982965150449277L,.313706645877887287337962201987L,.222381034453374470544355994426L,.101228536290376259152531354310L};
 long double sum=0,h=z/count;
 for(unsigned p=0;p<count;++p){auto mid=(p+.5L)*h;
  for(unsigned k=0;k<4;++k)for(int sign:{-1,1}){
   auto t=mid+sign*nodes[k]*h/2;
   auto value=1/std::sqrt(squared(t,om,w,wa));if(clock)value/=1+t;
   sum+=h/2*weights[k]*value;
  }
 }return sum;
}
}
int main(){try{
 cosmology::Policy policy;policy.integration={1e-14,1e-13,100000,30};
 for(double om:{0.,.3,1.})for(double w:{-2.,-1.,0.})for(double wa:{-2.,.7,2.}){
  auto model=cosmology::prepare_cpl({70,om,w,wa});check(model.status()==cosmology::Status::ok);
  check(model.model_id()=="P01/flat-cpl-radiation-free/v1");
  check(model.parameters().w0==w&&model.parameters().wa==wa);
  for(double z:{0.,1e-4,.5,5.}){
   cosmology::Query query{z,z,cosmology::Convention::geometric_same_redshift};
   auto batch=model.evaluate_batch(std::span(&query,1),policy);check(batch.slots.size()==1);
   const auto&s=batch.slots[0];check(s.status==cosmology::Status::ok);
   long double u=1+(long double)z, A=3*(1+(long double)w+wa)/u-3*(long double)wa/(u*u);
   auto Ap=-3*(1+(long double)w+wa)/(u*u)+6*(long double)wa/(u*u*u);
   auto D=std::pow(u,3*(1+(long double)w+wa))*std::exp(-3*(long double)wa*z/u);
   auto S=squared(z,om,w,wa),Sp=3*(long double)om*u*u+(1-om)*D*A;
   auto Spp=6*(long double)om*u+(1-om)*D*(A*A+Ap);
   near(s.expansion_E,std::sqrt(S),2e-13L,2e-12L);
   near(s.deceleration_q,-1+u*Sp/(2*S),2e-13L,2e-12L);
   near(s.jerk,1-u*Sp/S+u*u*Spp/(2*S),2e-13L,2e-12L);
   auto r16=panels(z,om,w,wa,16,false),r32=panels(z,om,w,wa,32,false);
   auto ref_ratio=std::abs(r32-r16)/(2e-15L+2e-10L*std::abs(r32));
   max_ref_ratio=std::max(max_ref_ratio,ref_ratio);check(ref_ratio<=.1L);
   near(s.radial_integral,r32,2e-15L,2e-10L);
   near(s.dimensionless_luminosity_shape,u*r32,u*2e-15L,2e-10L);
   auto clock=panels(z,om,w,wa,32,true);
   auto clock16=panels(z,om,w,wa,16,true);
   auto clock_ref_ratio=std::abs(clock-clock16)/(2e-15L+2e-10L*std::abs(clock));
   max_ref_ratio=std::max(max_ref_ratio,clock_ref_ratio);check(clock_ref_ratio<=.1L);
   // Hubble time from exact parsec definition, independent unit conversion.
   auto pc=648000.L*149597870700.L/acosl(-1.L),time=pc*1e6L/(70*1000.L);
   near(s.lookback_seconds/time,clock,2e-15L,2e-10L);
  }
 }
 // Constant-w all-dark exact power law, including p=0 and p=1 both sides.
 for(double w:{-1.,-1.+1e-9,-1.-1e-9,-1./3,-1./3+1e-9,-1./3-1e-9,0.}){
  auto g=cosmology::prepare_cpl({70,0,w,0});double z=1.3;
  cosmology::Query q{z,z,cosmology::Convention::geometric_same_redshift};
  auto s=g.evaluate_batch(std::span(&q,1),policy).slots[0];check(s.status==cosmology::Status::ok);
  auto p=1.5L*(1+(long double)w),L=log1pl(z),t=1-p;
  auto radial=t==0?L:expm1l(t*L)/t;
  near(s.radial_integral,radial,2e-15L,2e-10L);
 }
 auto legacy=cosmology::prepare({cosmology::Model::flat_cpl_late_v1,70,.3,0,-1,0});
 check(legacy.status()!=cosmology::Status::ok);
 for(auto p:{cosmology::CplParameters{70,-.1,-1,0}, {70,.3,-2.1,0},{70,.3,-1,2.1},{0,.3,-1,0},{70,.3,std::numeric_limits<double>::quiet_NaN(),0}})
  check(cosmology::prepare_cpl(p).status()!=cosmology::Status::ok);
 auto g=cosmology::prepare_cpl({70,.3,-.8,.7});
 cosmology::Query q{1,1,cosmology::Convention::geometric_same_redshift};
 auto tight=policy;tight.maximum_total_evaluations=1;
 check(g.evaluate_batch(std::span(&q,1),tight).slots[0].status==cosmology::Status::work_limit);
 q.z_expansion=6;check(g.evaluate_batch(std::span(&q,1),policy).slots[0].status==cosmology::Status::unsupported_domain);
 q={1,.9,cosmology::Convention::geometric_same_redshift};
 check(g.evaluate_batch(std::span(&q,1),policy).slots[0].status==cosmology::Status::incompatible_convention);
 std::printf("{\"suite\":\"background_cpl_peer\",\"checks\":%u,\"max_budget_fraction\":%.17Lg,\"max_reference_refinement_fraction\":%.17Lg,\"passed\":true}\n",checks,max_ratio,max_ref_ratio);
 return 0;
}catch(const std::exception&e){std::fprintf(stderr,"%s checks=%u\n",e.what(),checks);return 1;}}
