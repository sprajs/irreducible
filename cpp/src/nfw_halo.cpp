#include "irred/nfw_halo.hpp"
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
namespace irred::lensing {
namespace {
using W = long double;
using S = numerics::Status;
bool positive(double x) { return std::isfinite(x) && x > 0; }
bool valid(NFWHalo h, double r, double tol) {
  return positive(h.rho_s_msun_mpc3) && positive(h.radius_s_mpc) &&
         std::isfinite(r) && r >= 0 && positive(tol) &&
         std::fegetround() == FE_TONEAREST && std::numeric_limits<W>::digits >= 64;
}
numerics::ScalarResult output(W v, W scale, double tol) {
  double d = static_cast<double>(v);
  W e = 512 * std::numeric_limits<W>::epsilon() * scale +
        16 * std::numeric_limits<double>::epsilon() * std::abs(v) + std::abs(v-d);
  if (!std::isfinite(v) || !std::isfinite(d) || (v != 0 && !std::isnormal(d)))
    return {S::outside_domain};
  if (e > tol * std::abs(v) || !std::isfinite(static_cast<double>(e)))
    return {S::conditioning_budget_exceeded};
  return {S::ok,d,static_cast<double>(e),0};
}
// f = Sigma/(2 rho_s r_s), g = M_cyl/(4 pi rho_s r_s^3).
void shapes(W x, W &f, W &g) {
  W t=x-1;
  if (x < 0.01L) {
    W l=std::log(2/x), z=x*x;
    f=l-1+z*(1.5L*l-1.25L+z*(1.875L*l-47.L/32+z*(35.L/16*l-319.L/192)));
    g=z*(l/2-0.25L+z*(3.L/8*l-7.L/32+z*(5.L/16*l-37.L/192+z*(35.L/128*l-533.L/3072))));
  } else if (std::abs(t) < 0.001L) {
    f=1.L/3+t*(-2.L/5+t*(13.L/35+t*(-20.L/63+t*61.L/231)));
    // g'=x f; same local potential derivative identity.
    g=1-std::log(2.L)+t*(1.L/3+t*(-1.L/30+t*(-1.L/105+t*(17.L/1260+t*(-37.L/3465+t*61.L/1386)))));
  } else {
    W a;
    if (x<1) a=std::acosh(1/x)/std::sqrt((1-x)*(1+x));
    else a=std::acos(1/x)/std::sqrt((x-1)*(x+1));
    f=(1-a)/((x-1)*(x+1));
    g=std::log(x/2)+a;
  }
}
}
numerics::ScalarResult nfw_enclosed_mass(NFWHalo h, double r, double tol) noexcept {
  if (!valid(h,r,tol)) return {};
  if (r==0) return {S::ok,0,0,0};
  W x=W(r)/h.radius_s_mpc;
  if (x<1e-8L || x>1e8L) return {S::outside_domain};
  W m;
  if (x<0.001L) {
    // log(1+x)-x/(1+x), original alternating power series.
    m=0; W power=x*x;
    for (int n=2;n<=12;++n) { m+=(n%2? -1:1)*W(n-1)/n*power; power*=x; }
  } else m=std::log1p(x)-x/(1+x);
  W pref=4*std::numbers::pi_v<W>*h.rho_s_msun_mpc3*W(h.radius_s_mpc)*h.radius_s_mpc*h.radius_s_mpc;
  return output(pref*m,pref*std::abs(m),tol);
}
NFWProjection nfw_projection(NFWHalo h, double r, double tol) noexcept {
  NFWProjection out;
  if (!valid(h,r,tol)) return out;
  if (r==0) { out.status=S::singular; return out; }
  W x=W(r)/h.radius_s_mpc;
  if (x<1e-8L || x>1e8L) {out.status=S::outside_domain; return out;}
  W f,g; shapes(x,f,g);
  W pref=2*W(h.rho_s_msun_mpc3)*h.radius_s_mpc, mean=2*g/(x*x);
  // The local series truncation and special-function cancellation are included
  // in the scale diagnostic; this is not a certified libm enclosure.
  W scale=std::abs(f)+std::abs(mean);
  if (std::abs(x-1)<0.001L) scale+=std::pow(std::abs(x-1),5)/(512*std::numeric_limits<W>::epsilon());
  if (x<0.01L) scale+=4*std::pow(x,8)*(std::log(2/x)+1)/(512*std::numeric_limits<W>::epsilon());
  else if (x<0.5L) scale+=4*(std::abs(std::log(x/2))+1)/(x*x);
  if (std::abs(x-1)>=0.001L && x>=0.01L && x<2) scale+=1/std::abs(x-1);
  out.surface_density=output(pref*f,pref*scale,tol);
  out.mean_surface_density=output(pref*mean,pref*scale,tol);
  out.excess_surface_density=output(pref*(mean-f),pref*scale,tol);
  out.status=S::ok;
  for (auto v : {out.surface_density,out.mean_surface_density,out.excess_surface_density})
    if (v.status!=S::ok) {out.status=v.status; break;}
  return out;
}
NFWLens nfw_lens(NFWHalo h, double r, double critical, double distance, double tol) noexcept {
  NFWLens out;
  if (!positive(critical)||!positive(distance)) return out;
  auto p=nfw_projection(h,r,tol);
  if(p.status!=S::ok) {out.status=p.status; return out;}
  W k=W(p.surface_density.value)/critical, mean=W(p.mean_surface_density.value)/critical;
  W scale=(std::abs(k)+std::abs(mean));
  // Include inherited projection diagnostics in the admission scale.
  W inherited=(p.surface_density.error_estimate+p.mean_surface_density.error_estimate)/critical;
  auto project=[&](W v,W s) {
    auto result=output(v,s,tol);
    if (result.status==S::ok) {
      W error=result.error_estimate+inherited*s/scale;
      if(error>tol*std::abs(v)) return numerics::ScalarResult{S::conditioning_budget_exceeded};
      result.error_estimate=static_cast<double>(error);
    }
    return result;
  };
  out.convergence=project(k,scale);
  out.tangential_shear=project(mean-k,scale);
  out.deflection_radians=project(mean*W(r)/distance,scale*W(r)/distance);
  out.potential_radial_curvature=project(2*k-mean,2*scale);
  out.potential_tangential_curvature=project(mean,scale);
  out.status=S::ok;
  for(auto v:{out.convergence,out.tangential_shear,out.deflection_radians,out.potential_radial_curvature,out.potential_tangential_curvature})
    if(v.status!=S::ok) {out.status=v.status; break;}
  return out;
}
}
