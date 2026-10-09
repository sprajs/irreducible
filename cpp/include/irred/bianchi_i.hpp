#pragma once
#include "irred/numerics.hpp"
#include <array>
#include <span>
#include <string_view>
#include <vector>
namespace irred::cosmology {
// Diagonal Bianchi I, isotropic perfect-fluid dust/radiation/Lambda stress.
// a=1 and beta_i=0 at the supplied comoving observer/normalization anchor.
// H_anchor is in inverse seconds, not necessarily a present-day H0.
struct BianchiIModel {
  double anchor_hubble_per_second = 0;
  double matter_fraction = 0, radiation_fraction = 0, lambda_fraction = 0;
  // Canonical trace-free ownership: s_z=-(long double(s_x)+s_y).
  // s_i=dot(beta_i)/H_anchor at a=1; Omega_shear=sum(s_i^2)/6.
  double shear_x_over_anchor_hubble = 0, shear_y_over_anchor_hubble = 0;
};
struct BianchiIQuery {
  double scale_factor = 1;
  // Components in the observer's orthonormal principal-axis frame.
  // Must already be unit length; never normalized by the implementation.
  std::array<double,3> observed_direction{1,0,0};
};
struct BianchiIPolicy {
  double absolute_beta_tolerance = 1e-10, relative_tolerance = 1e-9;
  unsigned maximum_depth = 30;
  std::size_t maximum_callbacks = 1000000, maximum_rows = 4096;
  std::size_t maximum_native_bytes = 4 * 1024 * 1024;
};
struct BianchiIRow {
  BianchiIQuery query{};
  long double expansion_over_anchor_hubble = 0, mean_hubble_per_second = 0;
  std::array<long double,3> beta{}, directional_scale_factors{}, axis_hubble_per_second{};
  long double shear_fraction = 0;
  long double one_plus_directional_redshift = 0, directional_redshift = 0;
  // Empirical quadrature + arithmetic diagnostics, not proven error bounds.
  long double maximum_beta_estimate = 0, redshift_relative_estimate = 0;
  long double beta_trace_residual = 0;
};
struct BianchiITrajectory {
  numerics::Status status = numerics::Status::invalid_input;
  BianchiIModel initial_state{};
  std::array<long double,3> anchor_shear_over_hubble{};
  long double anchor_shear_fraction = 0;
  std::vector<BianchiIRow> rows;
  std::size_t callbacks = 0;
};
// Queries ordered from observer toward the past: nonincreasing a in [1e-4,1].
// Repeated epochs reuse the same integrated anisotropy, retaining every row and
// its supplied direction. Any refusal withholds all rows; no borrowed inputs.
BianchiITrajectory evolve_bianchi_i(const BianchiIModel &,
                                  std::span<const BianchiIQuery>, BianchiIPolicy = {});
inline constexpr std::string_view bianchi_i_model_id =
    "diagonal-Bianchi-I-GR-isotropic-perfect-fluid-dust-radiation-Lambda-shear/v1";
inline constexpr std::string_view bianchi_i_method_id =
    "shared-adaptive-Simpson-log-scale-scaled-shear-integral/v1";
} // namespace irred::cosmology
