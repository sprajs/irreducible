#pragma once
#include "irred/numerics.hpp"
#include <optional>
#include <span>
#include <string_view>
#include <vector>
namespace irred::cosmology {
// Flat GR: separately conserved ordinary dust, perfect-fluid radiation and
// a unified generalized Chaplygin barotropic fluid. All fractions at a=1.
struct ChaplyginModel {
  double anchor_hubble_km_s_mpc = 70;
  double ordinary_matter_fraction = 0, radiation_fraction = 0;
  double a_s = .7, alpha = 1;
  // Omega_cg is owned canonically as 1-(wide Omega_ordinary+Omega_radiation).
};
struct ChaplyginPolicy {
  double absolute_distance_tolerance_mpc = 1e-7, relative_tolerance = 1e-9;
  unsigned maximum_depth = 30;
  std::size_t maximum_callbacks = 1000000, maximum_rows = 4096;
  std::size_t maximum_native_bytes = 4 * 1024 * 1024;
};
struct ChaplyginDistances {
  long double comoving_mpc = 0, angular_diameter_mpc = 0, luminosity_mpc = 0;
  long double comoving_estimate_mpc = 0, angular_diameter_estimate_mpc = 0;
  long double luminosity_estimate_mpc = 0;
};
struct ChaplyginRow {
  double scale_factor = 1;
  long double expansion_over_anchor_hubble = 0, hubble_km_s_mpc = 0;
  long double fluid_density_over_anchor_fluid_density = 0;
  long double ordinary_matter_fraction = 0, radiation_fraction = 0, fluid_fraction = 0;
  long double fluid_equation_of_state = 0, fluid_one_plus_equation_of_state = 0;
  // Formal barotropic dp/de. This is not a collisionless signal speed or a
  // qualified rest-frame propagation law, especially on the vacuum endpoint.
  long double fluid_barotropic_slope = 0;
  long double total_equation_of_state = 0, deceleration = 0;
  long double background_relative_estimate = 0;
  std::optional<ChaplyginDistances> distances;
};
struct ChaplyginTrajectory {
  numerics::Status status = numerics::Status::invalid_input;
  ChaplyginModel initial_state{};
  long double anchor_fluid_fraction = 0;
  unsigned requested_outputs = 0;
  std::vector<ChaplyginRow> rows;
  std::size_t callbacks = 0;
};
inline constexpr unsigned chaplygin_distances = 1;
// Nonincreasing scales in [1e-4,1], duplicates retained/reused. Outputs are
// background only, or background plus all three flat distances. Atomic refusal
// withholds all rows; the batch owns its initial state and emitted values.
ChaplyginTrajectory evolve_chaplygin(const ChaplyginModel &,
                                   std::span<const double> scale_factors,
                                   unsigned requested_outputs = 0,
                                   ChaplyginPolicy = {});
inline constexpr std::string_view chaplygin_model_id =
    "flat-GR-generalized-Chaplygin-unified-barotropic-fluid-ordinary-dust-radiation/v1";
inline constexpr std::string_view chaplygin_method_id =
    "closed-log-sum-fluid-shared-adaptive-Simpson-flat-distance/v1";
} // namespace irred::cosmology
