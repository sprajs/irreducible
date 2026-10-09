#pragma once
#include "irred/numerics.hpp"
#include <span>
#include <string_view>
#include <vector>
namespace irred::cosmology {
// Flat GR homogeneous initial-value problem. H_i is at the supplied finite
// initial epoch, in inverse seconds; it is not a present-day H0.
struct DecayingMatterModel {
  double initial_scale_factor = 1;
  double initial_hubble_per_second = 0;
  double stable_matter_fraction = 0;
  double parent_fraction = 1;
  double daughter_radiation_fraction = 0;
  double lambda_fraction = 0;
  double decay_rate_over_initial_hubble = 0;
};
struct DecayingMatterPolicy {
  double absolute_tolerance = 1e-11;
  double relative_tolerance = 1e-9;
  std::size_t maximum_rhs_evaluations = 1000000;
  std::size_t maximum_rows = 4096;
  std::size_t maximum_native_bytes = 4 * 1024 * 1024;
};
struct DecayingMatterRow {
  double scale_factor = 0;
  long double elapsed_seconds = 0, elapsed_initial_hubble_time = 0;
  long double expansion_over_initial_hubble = 0, hubble_per_second = 0;
  long double stable_matter_fraction = 0, parent_fraction = 0;
  long double daughter_radiation_fraction = 0, lambda_fraction = 0;
  long double total_equation_of_state = 0, deceleration = 0;
  // Q=Gamma*rho_parent: Q/(H*rho_total), dimensionless and nonnegative.
  long double transfer_over_hubble_total_density = 0;
  // Densities in units rho_ci=3 H_i^2/(8 pi G), without supplying a G.
  long double stable_density_over_initial_critical = 0;
  long double parent_density_over_initial_critical = 0;
  long double daughter_density_over_initial_critical = 0;
  long double lambda_density_over_initial_critical = 0;
  long double total_density_over_initial_critical = 0;
  // Relative algebraic continuity residual (transfer cancels between species).
  long double continuity_residual = 0;
  // Accumulated local step-doubling estimates, not rigorous error bounds.
  long double elapsed_time_estimate = 0, comoving_daughter_estimate = 0;
  long double expansion_relative_estimate = 0;
};
struct DecayingMatterTrajectory {
  numerics::Status status = numerics::Status::invalid_input;
  DecayingMatterModel initial_state{};
  std::vector<DecayingMatterRow> rows;
  std::size_t rhs_evaluations = 0, accepted_steps = 0, rejected_steps = 0;
};
// Nondecreasing scales, duplicates retained, a_i<=a<=a_i exp(16).
// Any refusal withholds all rows. Owner retains the supplied initial state and
// outputs; no input spans/callbacks survive the call.
DecayingMatterTrajectory evolve_decaying_matter(
    const DecayingMatterModel &, std::span<const double> scale_factors,
    DecayingMatterPolicy = {});
inline constexpr std::string_view decaying_matter_model_id =
    "flat-GR-finite-anchor-pressureless-decay-massless-daughters/v1";
inline constexpr std::string_view decaying_matter_method_id =
    "positive-comoving-daughter-log-scale-RK4-step-doubling/v1";
} // namespace irred::cosmology
