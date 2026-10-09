#include "irred/plummer_sphere.hpp"
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
namespace irred::gravity {
namespace {
using W=long double;
using S=numerics::Status;
bool positive(double x) { return std::isfinite(x) && x>0; }
bool valid(PlummerSphere p,double radius,double tolerance) {
  return positive(p.total_mass_kg) && positive(p.scale_radius_metres) &&
         positive(p.gravitational_coupling_m3_kg_s2) && std::isfinite(radius) && radius>=0 &&
         positive(tolerance) && std::fegetround()==FE_TONEAREST &&
         std::numeric_limits<W>::digits>=64 && std::numeric_limits<W>::max_exponent>=16384;
}
numerics::ScalarResult project(W value,W scale,double tolerance) {
  double v=static_cast<double>(value);
  if (!std::isfinite(value) || !std::isfinite(v) ||
      (value!=0 && !std::isnormal(v))) return {S::outside_domain};
  W error=128*std::numeric_limits<W>::epsilon()*scale+
          std::abs(value-v)+8*std::numeric_limits<double>::epsilon()*std::abs(value);
  if(error>tolerance*std::abs(value) || !std::isfinite(static_cast<double>(error))) return {S::conditioning_budget_exceeded};
  return {S::ok,v,static_cast<double>(error),0};
}
}
PlummerEvaluation evaluate_plummer_sphere(PlummerSphere p,double radius,double tolerance) noexcept {
  PlummerEvaluation out;
  if(!valid(p,radius,tolerance)) return out;
  W r=radius,b=p.scale_radius_metres,s=std::hypot(r,b),q=r/s,t=b/s;
  W gm=W(p.gravitational_coupling_m3_kg_s2)*p.total_mass_kg;
  W mass=W(p.total_mass_kg)*q*q*q;
  W rho=3*W(p.total_mass_kg)*t*t/(4*std::numbers::pi_v<W>*s*s*s);
  W phi=-gm/s,acceleration=-gm*q/(s*s),speed=gm*q*q/s;
  W tangential=gm/(s*s*s),radial=tangential*(t*t-2*q*q);
  out.enclosed_mass_kg=project(mass,std::abs(mass),tolerance);
  out.density_kg_m3=project(rho,std::abs(rho),tolerance);
  out.potential_m2_s2=project(phi,std::abs(phi),tolerance);
  out.inward_radial_acceleration_m_s2=project(acceleration,std::abs(acceleration),tolerance);
  out.circular_speed_squared_m2_s2=project(speed,std::abs(speed),tolerance);
  out.potential_radial_curvature_s2=project(radial,tangential*(t*t+2*q*q),tolerance);
  out.potential_tangential_curvature_s2=project(tangential,tangential,tolerance);
  out.status=S::ok;
  for(auto v:{out.enclosed_mass_kg,out.density_kg_m3,out.potential_m2_s2,
              out.inward_radial_acceleration_m_s2,out.circular_speed_squared_m2_s2,
              out.potential_radial_curvature_s2,out.potential_tangential_curvature_s2})
    if(v.status!=S::ok) {out.status=v.status; break;}
  return out;
}
PlummerProjection project_plummer_sphere(PlummerSphere p,double radius,double tolerance) noexcept {
  PlummerProjection out;
  if(!valid(p,radius,tolerance)) return out;
  W r=radius,b=p.scale_radius_metres,s=std::hypot(r,b),q=r/s,t=b/s;
  W mass=W(p.total_mass_kg)*q*q,mean=W(p.total_mass_kg)/(std::numbers::pi_v<W>*s*s);
  W sigma=mean*t*t,excess=mean*q*q;
  out.projected_mass_kg=project(mass,mass,tolerance);
  out.surface_density_kg_m2=project(sigma,sigma,tolerance);
  out.mean_surface_density_kg_m2=project(mean,mean,tolerance);
  out.excess_surface_density_kg_m2=project(excess,excess,tolerance);
  out.status=S::ok;
  for(auto v:{out.projected_mass_kg,out.surface_density_kg_m2,
              out.mean_surface_density_kg_m2,out.excess_surface_density_kg_m2})
    if(v.status!=S::ok) {out.status=v.status; break;}
  return out;
}
PlummerLens plummer_thin_lens(PlummerSphere p,double radius,double critical,double distance,
                              double tolerance) noexcept {
  PlummerLens out;
  if(!valid(p,radius,tolerance) || !positive(critical) || !positive(distance)) return out;
  W r=radius,b=p.scale_radius_metres,s=std::hypot(r,b),q=r/s,t=b/s;
  // Independently supplied critical density includes the caller's geometry/G.
  W mean=W(p.total_mass_kg)/(std::numbers::pi_v<W>*W(critical)*s*s);
  W kappa=mean*t*t,shear=mean*q*q,alpha=mean*r/distance;
  W ratio=r/b;
  W potential=W(p.total_mass_kg)*std::log1p(ratio*ratio)/
              (2*std::numbers::pi_v<W>*W(critical)*distance*distance);
  // Factor the zero at r=b instead of subtracting squared radii.
  W radial=mean*((b-r)/s)*((b+r)/s);
  out.convergence=project(kappa,kappa,tolerance);
  out.tangential_shear=project(shear,shear,tolerance);
  out.deflection_radians=project(alpha,alpha,tolerance);
  out.potential_radians_squared=project(potential,potential,tolerance);
  out.potential_radial_curvature=project(radial,std::abs(radial),tolerance);
  out.potential_tangential_curvature=project(mean,mean,tolerance);
  out.status=S::ok;
  for(auto v:{out.convergence,out.tangential_shear,out.deflection_radians,
              out.potential_radians_squared,out.potential_radial_curvature,out.potential_tangential_curvature})
    if(v.status!=S::ok) {out.status=v.status; break;}
  return out;
}
}
