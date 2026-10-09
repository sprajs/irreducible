#include "irred/chaplygin.hpp"
#include "irred/quantities.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
using W=long double;
using S=numerics::Status;
constexpr W eps=std::numeric_limits<double>::epsilon();
struct State { W f=0,w=0,one_plus_w=0,e=0,matter=0,radiation=0,fluid=0,estimate=0; };
struct Context {
 const ChaplyginModel &m;
 W omega;
 mutable S failure=S::ok;
 bool state(W x,State &s) const {
  W logf=0;
  if(m.a_s==0) { logf=-3*x;s.w=0;s.one_plus_w=1; }
  else if(m.a_s==1) { logf=0;s.w=-1;s.one_plus_w=0; }
  else {
   const W first=std::log(W(m.a_s));
   const W second=std::log1p(-W(m.a_s))-3*(1+W(m.alpha))*x;
   const W high=std::max(first,second),low=std::min(first,second);
   const W logb=high+std::log1p(std::exp(low-high));
   logf=logb/(1+W(m.alpha));
   s.w=-std::exp(first-logb);s.one_plus_w=std::exp(second-logb);
  }
  // F(1)=1 is the defining normalization, not a correction to input closure.
  s.f=x==0?1:std::exp(logf);
  s.matter=m.ordinary_matter_fraction*std::exp(-3*x);
  s.radiation=m.radiation_fraction*std::exp(-4*x);
  s.fluid=omega*s.f;
  s.e=std::sqrt((s.matter+s.radiation)+s.fluid);
  s.estimate=128*std::numeric_limits<W>::epsilon()*(1+std::abs(logf)+std::abs(x));
  if(!std::isnormal(s.f) || !std::isnormal(s.fluid) || !std::isnormal(s.e) ||
     (m.ordinary_matter_fraction>0 && !std::isnormal(s.matter)) ||
     (m.radiation_fraction>0 && !std::isnormal(s.radiation)) ||
     (m.a_s>0 && !std::isnormal(s.w)) || (m.a_s<1 && !std::isnormal(s.one_plus_w))) {
   failure=S::outside_domain;return false;
  }
  return true;
 }
};
double distance_integrand(double x,const void *raw) {
 const auto &c=*static_cast<const Context *>(raw);State s;
 if(!c.state(x,s))return std::numeric_limits<double>::quiet_NaN();
 const double f=static_cast<double>(std::exp(-W(x))/s.e);
 if(!std::isnormal(f) || !(f>0)) { c.failure=S::outside_domain;return std::numeric_limits<double>::quiet_NaN(); }
 return f;
}
} // namespace
ChaplyginTrajectory evolve_chaplygin(const ChaplyginModel &m,
 std::span<const double> scales,unsigned requested,ChaplyginPolicy p) {
 ChaplyginTrajectory out;out.initial_state=m;out.requested_outputs=requested;
 const double inputs[]{m.anchor_hubble_km_s_mpc,m.ordinary_matter_fraction,m.radiation_fraction,
   m.a_s,m.alpha,p.absolute_distance_tolerance_mpc,p.relative_tolerance};
 for(double v:inputs)if(!std::isfinite(v)){out.status=S::nonfinite_input;return out;}
 if(std::fegetround()!=FE_TONEAREST || requested>chaplygin_distances ||
    !(p.absolute_distance_tolerance_mpc>0) || p.relative_tolerance<128*eps || p.maximum_depth>60) return out;
 if(m.anchor_hubble_km_s_mpc<1 || m.anchor_hubble_km_s_mpc>1000 ||
    m.ordinary_matter_fraction<0 || m.radiation_fraction<0 || m.a_s<0 || m.a_s>1 || m.alpha<0 || m.alpha>1) {
  out.status=S::outside_domain;return out;
 }
 out.anchor_fluid_fraction=1-(W(m.ordinary_matter_fraction)+m.radiation_fraction);
 if(!(out.anchor_fluid_fraction>0) || !std::isnormal(out.anchor_fluid_fraction)) {
  out.status=S::outside_domain;return out;
 }
 constexpr std::size_t overhead=sizeof(ChaplyginTrajectory)+32768;
 if(scales.size()>4096 || scales.size()>p.maximum_rows || p.maximum_native_bytes<overhead ||
    scales.size()>(p.maximum_native_bytes-overhead)/sizeof(ChaplyginRow)) {
  out.status=S::work_limit;return out;
 }
 double previous=1;
 for(double a:scales) {
  if(!std::isfinite(a)){out.status=S::nonfinite_input;return out;}
  if(a<1e-4 || a>previous){out.status=S::outside_domain;return out;}
  previous=a;
 }
 out.rows.reserve(scales.size());
 const bool distances=(requested&chaplygin_distances)!=0;
 const W unit=(W(irred::speed_of_light_m_per_s)/1000)/m.anchor_hubble_km_s_mpc;
 Context context{m,out.anchor_fluid_fraction};
 W chi=0,error=0;double prior_log=0;
 const double span=scales.empty()?0:-std::log(scales.back());
 const W minimum_a=scales.empty()?1:scales.back();
 for(double a:scales) {
  const double x=std::log(a);
  if(distances && x<prior_log) {
   if(p.maximum_callbacks-out.callbacks<3){out.status=S::work_limit;out.rows.clear();return out;}
   const numerics::IntegrationPolicy policy{
    static_cast<double>(p.absolute_distance_tolerance_mpc*minimum_a/unit/8*(prior_log-x)/std::max(1.0,span)),
    p.relative_tolerance/8,p.maximum_callbacks-out.callbacks,p.maximum_depth};
   const auto part=numerics::integrate(distance_integrand,&context,x,prior_log,policy);
   out.callbacks+=part.evaluations;
   if(context.failure!=S::ok || part.status!=S::ok) {
    out.status=context.failure!=S::ok?context.failure:part.status;out.rows.clear();return out;
   }
   if(!(part.value>0)){out.status=S::outside_domain;out.rows.clear();return out;}
   chi+=part.value;error+=part.error_estimate+64*eps*std::abs(part.value);
  }
  prior_log=x;
  State s;if(!context.state(std::log(W(a)),s)){out.status=context.failure;out.rows.clear();return out;}
  ChaplyginRow row;row.scale_factor=a;
  row.expansion_over_anchor_hubble=s.e;row.hubble_km_s_mpc=m.anchor_hubble_km_s_mpc*s.e;
  row.fluid_density_over_anchor_fluid_density=s.f;
  const W total=s.e*s.e;
  row.ordinary_matter_fraction=s.matter/total;row.radiation_fraction=s.radiation/total;row.fluid_fraction=s.fluid/total;
  row.fluid_equation_of_state=s.w;row.fluid_one_plus_equation_of_state=s.one_plus_w;
  row.fluid_barotropic_slope=-m.alpha*s.w;
  row.total_equation_of_state=row.fluid_fraction*s.w+row.radiation_fraction/3;
  row.deceleration=(1+3*row.total_equation_of_state)/2;row.background_relative_estimate=s.estimate;
  if(!std::isnormal(row.hubble_km_s_mpc) || !std::isnormal(row.fluid_fraction) ||
     (s.matter>0 && !std::isnormal(row.ordinary_matter_fraction)) ||
     (s.radiation>0 && !std::isnormal(row.radiation_fraction)) ||
     (m.alpha>0 && m.a_s>0 && !std::isnormal(row.fluid_barotropic_slope))) {
   out.status=S::outside_domain;out.rows.clear();return out;
  }
  if(s.estimate>p.relative_tolerance){out.status=S::conditioning_budget_exceeded;out.rows.clear();return out;}
  if(distances) {
   ChaplyginDistances d;
   d.comoving_mpc=unit*chi;d.angular_diameter_mpc=W(a)*d.comoving_mpc;d.luminosity_mpc=d.comoving_mpc/a;
   d.comoving_estimate_mpc=unit*error+s.estimate*d.comoving_mpc;
   d.angular_diameter_estimate_mpc=W(a)*d.comoving_estimate_mpc;
   d.luminosity_estimate_mpc=d.comoving_estimate_mpc/a;
   const W values[]{d.comoving_mpc,d.angular_diameter_mpc,d.luminosity_mpc};
   const W estimates[]{d.comoving_estimate_mpc,d.angular_diameter_estimate_mpc,d.luminosity_estimate_mpc};
   for(unsigned i=0;i<3;++i) {
    if(a<1 && (!(values[i]>0) || !std::isnormal(values[i]))) {out.status=S::outside_domain;out.rows.clear();return out;}
    if(!std::isfinite(estimates[i]) || estimates[i]>p.absolute_distance_tolerance_mpc+p.relative_tolerance*values[i]) {
     out.status=S::conditioning_budget_exceeded;out.rows.clear();return out;
    }
   }
   row.distances=d;
  }
  out.rows.push_back(row);
 }
 out.status=S::ok;return out;
}
} // namespace irred::cosmology
