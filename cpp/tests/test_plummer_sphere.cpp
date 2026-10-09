#include "irred/plummer_sphere.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
using namespace irred::gravity;
using S=irred::numerics::Status;
void close(double value,double reference,double relative=2e-9) {
  if(std::abs(value-reference)>relative*std::abs(reference)) {
    std::cerr<<std::setprecision(17)<<"comparison failure value="<<value<<" reference="<<reference<<" relative="<<relative<<'\n';
    std::abort();
  }
}
double density(double r) {return 3/(4*std::numbers::pi*std::pow(1+r*r,2.5));}
double radial_mass(double r,int n) {
  double step=std::atan(r)/n,sum=0;
  for(int i=0;i<=n;++i) {
    double u=i*step,t=std::tan(u),v=4*std::numbers::pi*density(t)*t*t/std::pow(std::cos(u),2);
    sum+=(i==0||i==n?1:i%2?4:2)*v;
  }
  return sum*step/3;
}
double shell_potential(double r,int n) {
  double lower=std::atan(r),step=(std::numbers::pi/2-lower)/n,sum=0;
  for(int i=0;i<=n;++i) {
    double u=lower+i*step,t=std::tan(u);
    double v=i==n?0:4*std::numbers::pi*density(t)*t/std::pow(std::cos(u),2);
    sum+=(i==0||i==n?1:i%2?4:2)*v;
  }
  return -(r==0?0:radial_mass(r,n)/r)-sum*step/3;
}
double surface(double r,int n) {
  double scale=std::hypot(r,1),step=(std::numbers::pi/2)/n,sum=0;
  for(int i=0;i<=n;++i) {
    double u=i*step,z=scale*std::tan(u);
    double v=i==n?0:2*density(std::hypot(r,z))*scale/std::pow(std::cos(u),2);
    sum+=(i==0||i==n?1:i%2?4:2)*v;
  }
  return sum*step/3;
}
double cylinder_mean(double r,int n) {
  double step=32.0/n,sum=0;
  for(int i=0;i<=n;++i) {
    double u=i*step,t=r*std::cosh(u),fraction=1/(std::cosh(u)*std::cosh(u)*(1+std::tanh(u)));
    double v=4*std::numbers::pi*density(t)*t*t*t*std::tanh(u)*fraction;
    sum+=(i==0||i==n?1:i%2?4:2)*v;
  }
  return (radial_mass(r,n)+sum*step/3)/(std::numbers::pi*r*r);
}
int main() {
  PlummerSphere model{1,1,1};
  for(double r:{1e-6,0.001,0.1,0.5,1.0,10.0,1000.0,1e6}) {
    auto e=evaluate_plummer_sphere(model,r); assert(e.status==S::ok);
    auto p=project_plummer_sphere(model,r); assert(p.status==S::ok);
    auto lens=plummer_thin_lens(model,r,10,100); assert(lens.status==S::ok);
    double m1=radial_mass(r,8192),m2=radial_mass(r,16384);
    double phi1=shell_potential(r,8192),phi2=shell_potential(r,16384);
    double s1=surface(r,8192),s2=surface(r,16384);
    double c1=cylinder_mean(r,8192),c2=cylinder_mean(r,16384);
    close(m1,m2,1e-10); close(phi1,phi2,1e-10);
    close(s1,s2,1e-10); close(c1,c2,1e-9);
    close(e.enclosed_mass_kg.value,m2); close(e.potential_m2_s2.value,phi2);
    close(p.surface_density_kg_m2.value,s2); close(p.mean_surface_density_kg_m2.value,c2,2e-8);
    close(-e.inward_radial_acceleration_m_s2.value,e.enclosed_mass_kg.value/(r*r),1e-14);
    close(e.circular_speed_squared_m2_s2.value,-r*e.inward_radial_acceleration_m_s2.value,1e-14);
    // Poisson residual is compared to inherited absolute output diagnostics,
    // since opposite curvatures in the far tail cancel to tiny density.
    double residual=e.potential_radial_curvature_s2.value+2*e.potential_tangential_curvature_s2.value-4*std::numbers::pi*e.density_kg_m3.value;
    double error=e.potential_radial_curvature_s2.error_estimate+2*e.potential_tangential_curvature_s2.error_estimate+4*std::numbers::pi*e.density_kg_m3.error_estimate;
    assert(std::abs(residual)<=error);
    close(lens.convergence.value,p.surface_density_kg_m2.value/10,1e-14);
    close(lens.tangential_shear.value,p.excess_surface_density_kg_m2.value/10,1e-14);
    close(lens.deflection_radians.value,p.mean_surface_density_kg_m2.value*r/1000,1e-14);
    double dr=r*1e-5;
    auto a=evaluate_plummer_sphere(model,r-dr),b=evaluate_plummer_sphere(model,r+dr);
    if(r<=1000) close((b.enclosed_mass_kg.value-a.enclosed_mass_kg.value)/(2*dr),4*std::numbers::pi*r*r*e.density_kg_m3.value,2e-5);
    if(r>=0.1) close((b.potential_m2_s2.value-a.potential_m2_s2.value)/(2*dr),-e.inward_radial_acceleration_m_s2.value,2e-7);
    close(-(b.inward_radial_acceleration_m_s2.value-a.inward_radial_acceleration_m_s2.value)/(2*dr),e.potential_radial_curvature_s2.value,2e-7);
    auto la=plummer_thin_lens(model,r-dr,10,100),lb=plummer_thin_lens(model,r+dr,10,100);
    if(r!=1) close((lb.deflection_radians.value-la.deflection_radians.value)*100/(2*dr),lens.potential_radial_curvature.value,2e-7);
    close((lb.potential_radians_squared.value-la.potential_radians_squared.value)*100/(2*dr),lens.deflection_radians.value,2e-7);
  }
  // Independent 90-/110-digit density-shell, LOS and cylindrical integrals
  // agreed to 1e-40 relative. Angular potential uses a Green-function integral
  // of the separately LOS-checked surface density instead of log1p.
  struct Fixture {double radius;std::array<double,17> value;};
  const Fixture fixtures[] = {
    {0.1, {0.00098518533684157340164878589479067,0.23286700428711438429707801611839,-0.99503719020998913566527375373857,-0.098518533684157340164878589479067,0.0098518533684157340164878589479067,0.95592240604429894417406948207411,0.98518533684157340164878589479067,0.0099009900990099009900990099009901,0.31203792391313662536787327393886,0.31515830315226799162155200667825,0.0031203792391313662536787327393886,0.031203792391313662536787327393886,0.00031203792391313662536787327393886,0.000031515830315226799162155200667825,0.000000015836443406814965898903574995115,0.030891754467400525911419454119947,0.031515830315226799162155200667825}},
    {1, {0.35355339059327376220042218105242,0.042202327319864347010399962078441,-0.70710678118654752440084436210485,-0.35355339059327376220042218105242,0.35355339059327376220042218105242,-0.17677669529663688110021109052621,0.35355339059327376220042218105242,0.5,0.079577471545947667884441881686257,0.15915494309189533576888376337251,0.079577471545947667884441881686257,0.0079577471545947667884441881686257,0.0079577471545947667884441881686257,0.00015915494309189533576888376337251,0.00000110317800076325796698228216059,0.0,0.015915494309189533576888376337251}},
    {10, {0.98518533684157340164878589479067,0.0000023286700428711438429707801611839,-0.099503719020998913566527375373857,-0.0098518533684157340164878589479067,0.098518533684157340164878589479067,-0.0019411077428858723458228553768648,0.00098518533684157340164878589479067,0.99009900990099009900990099009901,0.000031203792391313662536787327393886,0.0031515830315226799162155200667825,0.0031203792391313662536787327393886,0.0000031203792391313662536787327393886,0.00031203792391313662536787327393886,0.000031515830315226799162155200667825,0.0000073451924322010923746267778523877,-0.00030891754467400525911419454119947,0.00031515830315226799162155200667825}},
    {1e-12, {9.999999999999999999999985e-37,0.23873241463784300365332504822773,-0.9999999999999999999999995,-9.999999999999999999999985e-13,9.999999999999999999999985e-25,0.9999999999999999999999955,0.9999999999999999999999985,9.99999999999999999999999e-25,0.31830988618379067153776689012526,0.31830988618379067153776720843514,3.1830988618379067153776689012526e-25,0.031830988618379067153776689012526,3.1830988618379067153776689012526e-26,3.1830988618379067153776720843514e-16,1.5915494309189533576888368379504e-30,0.031830988618379067153776657181537,0.031830988618379067153776720843514}},
    {1e12, {0.9999999999999999999999985,2.3873241463784300365332504822773e-61,-9.999999999999999999999995e-13,-9.999999999999999999999985e-25,9.999999999999999999999985e-13,-1.999999999999999999999994e-36,9.999999999999999999999985e-37,0.999999999999999999999999,3.1830988618379067153776689012526e-49,3.1830988618379067153776720843514e-25,3.1830988618379067153776689012526e-25,3.1830988618379067153776689012526e-50,3.1830988618379067153776689012526e-26,3.1830988618379067153776720843514e-16,0.00008795227186553132890473449292026,-3.1830988618379067153776657181537e-26,3.1830988618379067153776720843514e-26}},
  };
  for(auto f:fixtures) {
    auto e=evaluate_plummer_sphere(model,f.radius);
    auto p=project_plummer_sphere(model,f.radius);
    auto l=plummer_thin_lens(model,f.radius,10,100);
    assert(e.status==S::ok && p.status==S::ok && l.status==S::ok);
    std::array values{e.enclosed_mass_kg,e.density_kg_m3,e.potential_m2_s2,
      e.inward_radial_acceleration_m_s2,e.circular_speed_squared_m2_s2,
      e.potential_radial_curvature_s2,e.potential_tangential_curvature_s2,
      p.projected_mass_kg,p.surface_density_kg_m2,p.mean_surface_density_kg_m2,p.excess_surface_density_kg_m2,
      l.convergence,l.tangential_shear,l.deflection_radians,l.potential_radians_squared,
      l.potential_radial_curvature,l.potential_tangential_curvature};
    for(unsigned i=0;i<values.size();++i) {
      if(f.value[i]==0) assert(values[i].value==0);
      else close(values[i].value,f.value[i],1e-14);
      assert(std::abs(values[i].value-f.value[i])<=values[i].error_estimate+
             2*std::numeric_limits<double>::epsilon()*std::abs(f.value[i]));
    }
  }
  auto scaled=evaluate_plummer_sphere({3,5,7},5);
  auto original=evaluate_plummer_sphere(model,1);
  assert(scaled.status==S::ok);
  close(scaled.enclosed_mass_kg.value,original.enclosed_mass_kg.value*3,1e-14);
  close(scaled.density_kg_m3.value,original.density_kg_m3.value*3/125,1e-14);
  close(scaled.potential_m2_s2.value,original.potential_m2_s2.value*21/5,1e-14);
  close(scaled.inward_radial_acceleration_m_s2.value,original.inward_radial_acceleration_m_s2.value*21/25,1e-14);
  close(scaled.circular_speed_squared_m2_s2.value,original.circular_speed_squared_m2_s2.value*21/5,1e-14);
  close(scaled.potential_radial_curvature_s2.value,original.potential_radial_curvature_s2.value*21/125,1e-14);
  close(scaled.potential_tangential_curvature_s2.value,original.potential_tangential_curvature_s2.value*21/125,1e-14);
  auto sp=project_plummer_sphere({3,5,7},5),op=project_plummer_sphere(model,1);
  assert(sp.status==S::ok);
  close(sp.projected_mass_kg.value,op.projected_mass_kg.value*3,1e-14);
  close(sp.surface_density_kg_m2.value,op.surface_density_kg_m2.value*3/25,1e-14);
  close(sp.mean_surface_density_kg_m2.value,op.mean_surface_density_kg_m2.value*3/25,1e-14);
  close(sp.excess_surface_density_kg_m2.value,op.excess_surface_density_kg_m2.value*3/25,1e-14);
  auto sl=plummer_thin_lens({3,5,7},5,1.2,500),ol=plummer_thin_lens(model,1,10,100);
  assert(sl.status==S::ok);
  close(sl.convergence.value,ol.convergence.value,1e-14);
  close(sl.tangential_shear.value,ol.tangential_shear.value,1e-14);
  close(sl.deflection_radians.value,ol.deflection_radians.value,1e-14);
  close(sl.potential_radians_squared.value,ol.potential_radians_squared.value,1e-14);
  close(sl.potential_tangential_curvature.value,ol.potential_tangential_curvature.value,1e-14);
  auto center=evaluate_plummer_sphere(model,0); assert(center.status==S::ok);
  assert(center.enclosed_mass_kg.value==0 && center.inward_radial_acceleration_m_s2.value==0 && center.circular_speed_squared_m2_s2.value==0);
  close(center.density_kg_m3.value,3/(4*std::numbers::pi),1e-14);
  assert(center.potential_m2_s2.value==-1 && center.potential_radial_curvature_s2.value==1 && center.potential_tangential_curvature_s2.value==1);
  auto p0=project_plummer_sphere(model,0); assert(p0.status==S::ok);
  assert(p0.projected_mass_kg.value==0 && p0.excess_surface_density_kg_m2.value==0);
  close(p0.surface_density_kg_m2.value,1/std::numbers::pi,1e-14);
  assert(p0.surface_density_kg_m2.value==p0.mean_surface_density_kg_m2.value);
  auto l0=plummer_thin_lens(model,0,10,100); assert(l0.status==S::ok);
  assert(l0.deflection_radians.value==0 && l0.potential_radians_squared.value==0 && l0.tangential_shear.value==0);
  assert(l0.potential_radial_curvature.value==l0.potential_tangential_curvature.value);
  auto l1=plummer_thin_lens(model,1,10,100); assert(l1.status==S::ok && l1.potential_radial_curvature.value==0);
  for(double r:{std::nextafter(1.0,0.0),std::nextafter(1.0,2.0)}) {
    auto lens=plummer_thin_lens(model,r,10,100); assert(lens.status==S::ok);
    assert((r<1)==(lens.potential_radial_curvature.value>0));
  }
  auto near=evaluate_plummer_sphere(model,1e-12),far=evaluate_plummer_sphere(model,1e12);
  assert(near.status==S::ok && far.status==S::ok);
  close(near.enclosed_mass_kg.value/1e-36,1,1e-14);
  close(near.inward_radial_acceleration_m_s2.value/1e-12,-1,1e-14);
  close(far.enclosed_mass_kg.value,1,1e-14);
  close(far.potential_m2_s2.value*1e12,-1,1e-14);
  close(far.inward_radial_acceleration_m_s2.value*1e24,-1,1e-14);
  auto near_projection=project_plummer_sphere(model,1e-12);
  close(near_projection.excess_surface_density_kg_m2.value/1e-24,1/std::numbers::pi,1e-14);
  auto tidal_zero=evaluate_plummer_sphere(model,std::sqrt(0.5));
  assert(tidal_zero.status==S::conditioning_budget_exceeded);
  assert(tidal_zero.potential_radial_curvature_s2.status==S::conditioning_budget_exceeded);
  assert(tidal_zero.potential_tangential_curvature_s2.status==S::ok);
  assert(evaluate_plummer_sphere({0,1,1},1).status==S::invalid_input);
  assert(evaluate_plummer_sphere({1,-1,1},1).status==S::invalid_input);
  assert(evaluate_plummer_sphere({1,1,0},1).status==S::invalid_input);
  assert(evaluate_plummer_sphere(model,-1).status==S::invalid_input);
  assert(evaluate_plummer_sphere(model,std::numeric_limits<double>::infinity()).status==S::invalid_input);
  assert(evaluate_plummer_sphere(model,1,1e-30).status==S::conditioning_budget_exceeded);
  assert(project_plummer_sphere(model,1,1e-30).status==S::conditioning_budget_exceeded);
  assert(plummer_thin_lens(model,1,0,100).status==S::invalid_input);
  assert(plummer_thin_lens(model,1,10,0).status==S::invalid_input);
  assert(evaluate_plummer_sphere({1e308,1e-308,1e308},1e-308).status==S::outside_domain);
  assert(project_plummer_sphere({1e308,1e-308,1},0).status==S::outside_domain);
  assert(plummer_thin_lens({1e308,1e-308,1},1e-308,1e-308,1e-308).status==S::outside_domain);
  assert(evaluate_plummer_sphere({1e-308,1e308,1e-308},1).status==S::outside_domain);
  assert(project_plummer_sphere({1e-308,1e308,1},1).status==S::outside_domain);
  assert(plummer_thin_lens({1e-308,1e308,1},1,1e308,1e308).status==S::outside_domain);
  std::fesetround(FE_DOWNWARD);
  assert(evaluate_plummer_sphere(model,1).status==S::invalid_input);
  assert(project_plummer_sphere(model,1).status==S::invalid_input);
  assert(plummer_thin_lens(model,1,10,100).status==S::invalid_input);
  std::fesetround(FE_TONEAREST);
}
