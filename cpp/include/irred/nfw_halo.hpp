#pragma once
#include "irred/numerics.hpp"
#include <string_view>
namespace irred::lensing {
inline constexpr std::string_view nfw_halo_id = "spherical-untruncated-NFW-supplied-critical-density/v1";
struct NFWHalo { double rho_s_msun_mpc3, radius_s_mpc; };
struct NFWProjection {
  numerics::Status status = numerics::Status::invalid_input;
  numerics::ScalarResult surface_density, mean_surface_density, excess_surface_density;
};
struct NFWLens {
  numerics::Status status = numerics::Status::invalid_input;
  numerics::ScalarResult convergence, tangential_shear, deflection_radians;
  // Eigenvalues of the angular potential Hessian, not the lens Jacobian.
  numerics::ScalarResult potential_radial_curvature, potential_tangential_curvature;
};
// Physical radii in Mpc; mass in Msun; surface densities in Msun/Mpc^2.
// Projection admits 1e-8 <= R/r_s <= 1e8. The origin is singular.
numerics::ScalarResult nfw_enclosed_mass(NFWHalo, double radius_mpc,
                                        double relative_tolerance = 1e-10) noexcept;
NFWProjection nfw_projection(NFWHalo, double radius_mpc,
                             double relative_tolerance = 1e-10) noexcept;
// Caller owns the thin-lens geometry and critical-density physical convention.
// Deflection = mean convergence * R / D_l. No distance model is inferred.
NFWLens nfw_lens(NFWHalo, double radius_mpc, double sigma_critical_msun_mpc2,
                 double lens_distance_mpc, double relative_tolerance = 1e-10) noexcept;
}
