#include "irred/hernquist_sphere.hpp"
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
namespace irred::gravity {
namespace {
using W=long double;
using S=numerics::Status;
bool positive(double v) { return std::isfinite(v) && v>0; }
numerics::ScalarResult project(W value, double tolerance) {
  double v=static_cast<double>(value);
  if (!std::isfinite(value) || !std::isfinite(v) ||
      (value!=0 && !std::isnormal(v))) return {S::outside_domain};
  // Straight products and positive sums; no cancellation or iterative solver.
  W error=64*std::numeric_limits<W>::epsilon()*std::abs(value)+
          std::abs(value-v)+8*std::numeric_limits<double>::epsilon()*std::abs(value);
  if (error>tolerance*std::abs(value)) return {S::conditioning_budget_exceeded};
  return {S::ok,v,static_cast<double>(error),0};
}
}
HernquistEvaluation evaluate_hernquist_sphere(HernquistSphere model, double r,
                                             double tolerance) noexcept {
  HernquistEvaluation out;
  if (!positive(model.total_mass_kg) || !positive(model.scale_radius_metres) ||
      !positive(model.gravitational_coupling_m3_kg_s2) || !std::isfinite(r) || r<0 ||
      !positive(tolerance) || std::fegetround()!=FE_TONEAREST ||
      std::numeric_limits<W>::digits<64 || std::numeric_limits<W>::max_exponent<16384)
    return out;
  W mass=model.total_mass_kg, a=model.scale_radius_metres, radius=r;
  W denominator=radius+a, gm=W(model.gravitational_coupling_m3_kg_s2)*mass;
  W q=radius/denominator;
  out.enclosed_mass_kg=project(mass*q*q,tolerance);
  out.potential_m2_s2=project(-gm/denominator,tolerance);
  out.circular_speed_squared_m2_s2=project(gm/denominator*q,tolerance);
  if (r==0) {
    out.density_kg_m3.status=S::singular;
    out.inward_radial_acceleration_m_s2.status=S::singular;
    out.potential_radial_curvature_s2.status=S::singular;
    out.potential_tangential_curvature_s2.status=S::singular;
  } else {
    W inverse_squared=1/(denominator*denominator);
    out.density_kg_m3=project(mass*(a/denominator)*inverse_squared/
                             (2*std::numbers::pi_v<W>*radius),tolerance);
    out.inward_radial_acceleration_m_s2=project(-gm*inverse_squared,tolerance);
    out.potential_radial_curvature_s2=project(-2*gm*inverse_squared/denominator,tolerance);
    out.potential_tangential_curvature_s2=project(gm*inverse_squared/radius,tolerance);
  }
  out.status=S::ok;
  for(auto value:{out.enclosed_mass_kg,out.density_kg_m3,out.potential_m2_s2,
                 out.inward_radial_acceleration_m_s2,out.circular_speed_squared_m2_s2,
                 out.potential_radial_curvature_s2,out.potential_tangential_curvature_s2})
    if(value.status!=S::ok) {out.status=value.status; break;}
  return out;
}
}
