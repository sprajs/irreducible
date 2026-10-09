#include "irred/nfw_halo.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <cfenv>
#include <limits>
#include <numbers>
using namespace irred::lensing;
using S=irred::numerics::Status;
void close(double got,double want,double tol=2e-9) { assert(std::abs(got-want)<=tol*std::abs(want)); }
// Original independent Simpson quadratures in logarithmic line-of-sight
// coordinates. Sigma=2∫rho(sqrt(R²+z²))dz, z=R sinh(u).
// M_cyl=4pi∫rho(r)r² fraction_inside_cylinder(r)dr, r=exp(t).
double sigma(double x,int n) {
  const double step=32.0/n; double sum=0;
  for(int i=0;i<=n;++i) {double u=i*step,r=x*std::cosh(u);
    double v=1/std::pow(1+r,2); sum+=(i==0||i==n?1:i%2?4:2)*v;}
  return 2*sum*step/3;
}
double cylinder(double x,int n) {
  // Split at the square-root cusp using r=x cosh(u) for the exterior.
  const double step=32.0/n; double sum=0;
  for(int i=0;i<=n;++i) {double u=i*step,r=x*std::cosh(u),t=std::tanh(u);
    double fraction=1/(std::cosh(u)*std::cosh(u)*(1+t));
    double v=r*r*t*fraction/std::pow(1+r,2);
    sum+=(i==0||i==n?1:i%2?4:2)*v;}
  double interior=std::log1p(x)-x/(1+x);
  return 4*(interior+sum*step/3)/(x*x);
}
int main() {
  NFWHalo h{1,1};
  for(double x:{0.00001,0.009999999,0.01,0.010000001,0.2,0.998999999,0.999,0.999000001,0.999999999,1.0,1.000000001,1.000999999,1.001,1.001000001,2.0,10.0,1000.0}) {
    auto p=nfw_projection(h,x); assert(p.status==S::ok);
    close(p.surface_density.value,sigma(x,16384));
    close(p.mean_surface_density.value,cylinder(x,16384),2e-8);
    close(sigma(x,8192),sigma(x,16384),1e-10);
    close(cylinder(x,8192),cylinder(x,16384),1e-9);
    close(p.surface_density.value,sigma(x,8192));
    close(p.mean_surface_density.value,cylinder(x,8192),2e-8);
    auto lens=nfw_lens(h,x,10,100); assert(lens.status==S::ok);
    close(lens.deflection_radians.value,p.mean_surface_density.value*x/1000);
    close(lens.tangential_shear.value,p.excess_surface_density.value/10);
    double dx=x*1e-5;
    auto a=nfw_lens(h,x-dx,10,100),b=nfw_lens(h,x+dx,10,100);
    close((b.deflection_radians.value-a.deflection_radians.value)*100/(2*dx),lens.potential_radial_curvature.value,2e-7);
    auto m=nfw_enclosed_mass(h,x); assert(m.status==S::ok);
    auto ma=nfw_enclosed_mass(h,x-dx),mb=nfw_enclosed_mass(h,x+dx);
    close((mb.value-ma.value)/(2*dx),4*std::numbers::pi*x/std::pow(1+x,2),2e-7);
  }
  // 60-digit mpmath direct LOS integrals, not the production closed forms.
  close(nfw_projection(h,1).surface_density.value,0.6666666666666666666666666667,1e-14);
  close(nfw_projection(h,1).mean_surface_density.value,1.2274112777602187623310715142,1e-14);
  close(nfw_projection(h,0.1).surface_density.value,4.057176060312765,1e-14);
  struct Fixture { double x,sigma,mean; };
  for(auto v : {Fixture{1e-8,36.227655849024627047,37.227655849024624343},
                Fixture{0.999999999,0.6666666674666666674,1.2274112788817079854},
                Fixture{1.000000001,0.6666666658666666674,1.2274112766387295410},
                Fixture{10,0.0172160855305560105,0.0702896671468631142},
                Fixture{1e8,1.9999999685840738641e-16,7.0910134316401533323e-15}}) {
    auto p=nfw_projection(h,v.x); assert(p.status==S::ok);
    close(p.surface_density.value,v.sigma,1e-13);
    close(p.mean_surface_density.value,v.mean,1e-13);
  }
  auto tiny=nfw_projection(h,1e-8); assert(tiny.status==S::ok);
  close(tiny.excess_surface_density.value,1.0,1e-13);
  auto large=nfw_projection(h,1e8); assert(large.status==S::ok);
  close(large.surface_density.value*1e16,2.0,2e-8);
  assert(nfw_enclosed_mass(h,0).status==S::ok);
  assert(nfw_projection(h,0).status==S::singular);
  assert(nfw_projection(h,1e-9).status==S::outside_domain);
  assert(nfw_projection({-1,1},1).status==S::invalid_input);
  assert(nfw_projection(h,std::numeric_limits<double>::infinity()).status==S::invalid_input);
  assert(nfw_projection({1e308,1e308},1e308).status==S::outside_domain);
  assert(nfw_projection(h,1,1e-30).status==S::conditioning_budget_exceeded);
  assert(nfw_lens(h,1,0,1).status==S::invalid_input);
  std::fesetround(FE_DOWNWARD); assert(nfw_projection(h,1).status==S::invalid_input); std::fesetround(FE_TONEAREST);
}
