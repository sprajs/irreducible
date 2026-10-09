#pragma once
#include "irred/numerics.hpp"
#include <string_view>
namespace irred::gravity {
inline constexpr std::string_view plummer_sphere_id =
    "finite-core-Newtonian-Plummer-SI-supplied-G-and-lens-geometry/v1";
struct PlummerSphere {
  double total_mass_kg, scale_radius_metres, gravitational_coupling_m3_kg_s2;
};
struct PlummerEvaluation {
  numerics::Status status = numerics::Status::invalid_input;
  numerics::ScalarResult enclosed_mass_kg, density_kg_m3, potential_m2_s2;
  // Signed outward radial component: negative means inward, zero at origin.
  numerics::ScalarResult inward_radial_acceleration_m_s2,
                         circular_speed_squared_m2_s2;
  // Cartesian potential Hessian eigenvalues, in s^-2.
  numerics::ScalarResult potential_radial_curvature_s2,
                         potential_tangential_curvature_s2;
};
struct PlummerProjection {
  numerics::Status status = numerics::Status::invalid_input;
  numerics::ScalarResult projected_mass_kg, surface_density_kg_m2,
                         mean_surface_density_kg_m2, excess_surface_density_kg_m2;
};
struct PlummerLens {
  numerics::Status status = numerics::Status::invalid_input;
  numerics::ScalarResult convergence, tangential_shear, deflection_radians;
  // Angular lens potential psi(0)=0, in radians^2, and its Hessian eigenvalues.
  numerics::ScalarResult potential_radians_squared,
                         potential_radial_curvature, potential_tangential_curvature;
};
// All radii in metres. Phi(infinity)=0. The origin is nonsingular.
PlummerEvaluation evaluate_plummer_sphere(PlummerSphere, double radius_metres,
                                          double relative_tolerance=1e-12) noexcept;
PlummerProjection project_plummer_sphere(PlummerSphere, double radius_metres,
                                         double relative_tolerance=1e-12) noexcept;
// Geometry, validity of thin-lens approximation and consistency of supplied
// critical density with G are caller-owned. No FLRW model is inferred.
PlummerLens plummer_thin_lens(PlummerSphere, double radius_metres,
                              double sigma_critical_kg_m2, double lens_distance_metres,
                              double relative_tolerance=1e-12) noexcept;
}
