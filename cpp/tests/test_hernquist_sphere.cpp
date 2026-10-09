#include "irred/hernquist_sphere.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <numbers>
using namespace irred::gravity;
using S=irred::numerics::Status;
void close(double value,double reference,double relative=2e-9) {
  if (std::abs(value-reference)>relative*std::abs(reference)) {
    std::cerr<<std::setprecision(17)<<"comparison failure value="<<value<<" reference="<<reference<<" relative="<<relative<<"\n";
    std::abort();
  }
}
double density(double r) {
  // Direct source density, separate from the production scaled algebra.
  return 1/(2*std::numbers::pi*r*std::pow(r+1,3));
}
double radial_mass(double r,int n) {
  double step=std::log1p(r)/n,sum=0;
  for(int i=0;i<=n;++i) {
    double u=i*step,t=std::expm1(u);
    double v=i==0?0:4*std::numbers::pi*density(t)*t*t*std::exp(u);
    sum+=(i==0||i==n?1:i%2?4:2)*v;
  }
  return sum*step/3;
}
double shell_potential(double r,int n) {
  // Newton's interior/exterior spherical-shell potential, density integrated
  // directly. Exterior radial map t=expm1(u), upper u=50 preserves the tail.
  double lower=std::log1p(r),step=(50-lower)/n,sum=0;
  for(int i=0;i<=n;++i) {
    double u=lower+i*step,t=std::expm1(u);
    double v=4*std::numbers::pi*density(t)*t*std::exp(u);
    sum+=(i==0||i==n?1:i%2?4:2)*v;
  }
  return -radial_mass(r,n)/r-sum*step/3;
}
int main() {
  HernquistSphere model{1,1,1};
  for(double r:{1e-6,0.001,0.1,1.0,10.0,1000.0,1e6}) {
    auto e=evaluate_hernquist_sphere(model,r); assert(e.status==S::ok);
    double m1=radial_mass(r,8192),m2=radial_mass(r,16384);
    double p1=shell_potential(r,8192),p2=shell_potential(r,16384);
    close(m1,m2,1e-10); close(p1,p2,1e-9);
    close(e.enclosed_mass_kg.value,m2);
    close(e.potential_m2_s2.value,p2,2e-8);
    close(-e.inward_radial_acceleration_m_s2.value,e.enclosed_mass_kg.value/(r*r),1e-14);
    close(e.circular_speed_squared_m2_s2.value,-r*e.inward_radial_acceleration_m_s2.value,1e-14);
    close(e.potential_radial_curvature_s2.value+2*e.potential_tangential_curvature_s2.value,
          4*std::numbers::pi*e.density_kg_m3.value,2e-9);
    // At tiny r a relative 1e-5 step makes binary64 potential differences
    // roundoff dominated. A wider centered step remains far below a and
    // retains the same derivative acceptance budgets.
    double dr=r*(r<0.01?0.01:1e-5);
    auto a=evaluate_hernquist_sphere(model,r-dr),b=evaluate_hernquist_sphere(model,r+dr);
    close((b.enclosed_mass_kg.value-a.enclosed_mass_kg.value)/(2*dr),
          4*std::numbers::pi*r*r*e.density_kg_m3.value,2e-5);
    close((b.potential_m2_s2.value-a.potential_m2_s2.value)/(2*dr),
          -e.inward_radial_acceleration_m_s2.value,2e-7);
    close(-(b.inward_radial_acceleration_m_s2.value-a.inward_radial_acceleration_m_s2.value)/(2*dr),
          e.potential_radial_curvature_s2.value,2e-7);
  }
  // Independent 70-/90-digit density-shell mpmath quadratures agreed to
  // 1e-50 relative; acceleration/tides then use shell gravity and Poisson.
  struct Fixture { double radius; std::array<double,7> value; };
  const Fixture fixtures[] = {
    {0.1, {0.008264462809917355371900826446281,1.1957546438158928307203889058791,-0.90909090909090909090909090909091,-0.8264462809917355371900826446281,0.08264462809917355371900826446281,-1.5026296018031555221637866265965,8.264462809917355371900826446281}},
    {1, {0.25,0.019894367886486916971110470421564,-0.5,-0.25,0.25,-0.25,0.25}},
    {10, {0.8264462809917355371900826446281,0.000011957546438158928307203889058791,-0.090909090909090909090909090909091,-0.008264462809917355371900826446281,0.08264462809917355371900826446281,-0.0015026296018031555221637866265965,0.0008264462809917355371900826446281}},
    {1e-12, {9.99999999998000000000003e-25,159154943091.41787093960903229487,-0.999999999999000000000001,-0.999999999998000000000003,9.99999999998000000000003e-13,-1.999999999994000000000012,999999999998.000000000003}},
    {1e12, {0.999999999998000000000003,1.5915494309141787093960903229487e-49,-9.99999999999000000000001e-13,-9.99999999998000000000003e-25,9.99999999998000000000003e-13,-1.999999999994000000000012e-36,9.99999999998000000000003e-37}},
  };
  for(auto fixture:fixtures) {
    auto e=evaluate_hernquist_sphere(model,fixture.radius); assert(e.status==S::ok);
    std::array values{e.enclosed_mass_kg,e.density_kg_m3,e.potential_m2_s2,
                      e.inward_radial_acceleration_m_s2,e.circular_speed_squared_m2_s2,
                      e.potential_radial_curvature_s2,e.potential_tangential_curvature_s2};
    for(unsigned i=0;i<values.size();++i) {
      close(values[i].value,fixture.value[i],1e-14);
      assert(std::abs(values[i].value-fixture.value[i])<=values[i].error_estimate+
             2*std::numeric_limits<double>::epsilon()*std::abs(fixture.value[i]));
    }
  }
  auto scaled=evaluate_hernquist_sphere({3,5,7},5);
  auto original=evaluate_hernquist_sphere(model,1);
  assert(scaled.status==S::ok);
  close(scaled.enclosed_mass_kg.value,original.enclosed_mass_kg.value*3,1e-14);
  close(scaled.density_kg_m3.value,original.density_kg_m3.value*3/125,1e-14);
  close(scaled.potential_m2_s2.value,original.potential_m2_s2.value*21/5,1e-14);
  close(scaled.inward_radial_acceleration_m_s2.value,original.inward_radial_acceleration_m_s2.value*21/25,1e-14);
  close(scaled.circular_speed_squared_m2_s2.value,original.circular_speed_squared_m2_s2.value*21/5,1e-14);
  close(scaled.potential_radial_curvature_s2.value,original.potential_radial_curvature_s2.value*21/125,1e-14);
  close(scaled.potential_tangential_curvature_s2.value,original.potential_tangential_curvature_s2.value*21/125,1e-14);
  // Exact rational identities at r=a, and the total-mass/Kepler tail.
  auto one=evaluate_hernquist_sphere(model,1);
  close(one.enclosed_mass_kg.value,0.25,1e-14);
  close(one.potential_m2_s2.value,-0.5,1e-14);
  close(one.inward_radial_acceleration_m_s2.value,-0.25,1e-14);
  close(one.circular_speed_squared_m2_s2.value,0.25,1e-14);
  close(one.potential_radial_curvature_s2.value,-0.25,1e-14);
  close(one.potential_tangential_curvature_s2.value,0.25,1e-14);
  auto near=evaluate_hernquist_sphere(model,1e-12);
  assert(near.status==S::ok);
  close(near.enclosed_mass_kg.value/1e-24,1,3e-12);
  close(near.inward_radial_acceleration_m_s2.value,-1,3e-12);
  auto far=evaluate_hernquist_sphere(model,1e12); assert(far.status==S::ok);
  close(far.enclosed_mass_kg.value,1,3e-12);
  close(far.potential_m2_s2.value*1e12,-1,2e-12);
  close(far.inward_radial_acceleration_m_s2.value*1e24,-1,3e-12);
  auto origin=evaluate_hernquist_sphere(model,0);
  assert(origin.status==S::singular);
  assert(origin.enclosed_mass_kg.status==S::ok && origin.enclosed_mass_kg.value==0);
  assert(origin.potential_m2_s2.status==S::ok && origin.potential_m2_s2.value==-1);
  assert(origin.circular_speed_squared_m2_s2.status==S::ok && origin.circular_speed_squared_m2_s2.value==0);
  assert(evaluate_hernquist_sphere({0,1,1},1).status==S::invalid_input);
  assert(evaluate_hernquist_sphere({1,0,1},1).status==S::invalid_input);
  assert(evaluate_hernquist_sphere({1,1,-1},1).status==S::invalid_input);
  assert(evaluate_hernquist_sphere(model,-1).status==S::invalid_input);
  assert(evaluate_hernquist_sphere(model,std::numeric_limits<double>::quiet_NaN()).status==S::invalid_input);
  assert(evaluate_hernquist_sphere(model,1,1e-30).status==S::conditioning_budget_exceeded);
  assert(evaluate_hernquist_sphere({1e308,1e-308,1e308},1e-308).status==S::outside_domain);
  assert(evaluate_hernquist_sphere({1e-308,1e308,1e-308},1).status==S::outside_domain);
  std::fesetround(FE_UPWARD);
  assert(evaluate_hernquist_sphere(model,1).status==S::invalid_input);
  std::fesetround(FE_TONEAREST);
}
