#pragma once
#include "irred/numerics.hpp"
#include <string_view>
namespace irred::gravity {
inline constexpr std::string_view hernquist_sphere_id =
    "finite-mass-Newtonian-Hernquist-supplied-G-SI/v1";
// SI throughout. G is supplied, not inferred or assigned a fixed asset value.
struct HernquistSphere {
  double total_mass_kg, scale_radius_metres, gravitational_coupling_m3_kg_s2;
};
struct HernquistEvaluation {
  numerics::Status status = numerics::Status::invalid_input;
  numerics::ScalarResult enclosed_mass_kg, density_kg_m3;
  // Acceleration is the signed outward radial component: negative is inward.
  numerics::ScalarResult potential_m2_s2, inward_radial_acceleration_m_s2;
  numerics::ScalarResult circular_speed_squared_m2_s2;
  // Eigenvalues of the Cartesian Hessian of the Newtonian potential.
  numerics::ScalarResult potential_radial_curvature_s2,
                         potential_tangential_curvature_s2;
};
// Phi(infinity)=0. Density, force direction and tidal curvatures are singular
// at the origin; the other three scalars retain their finite limiting values.
// Aggregate status reports any refused output; each scalar has its own status.
HernquistEvaluation evaluate_hernquist_sphere(
    HernquistSphere, double radius_metres,
    double relative_tolerance = 1e-12) noexcept;
}
