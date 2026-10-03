#pragma once
#include "irred/continuous_cmb_projection.hpp"
#include "irred/thermal_neutrino.hpp"
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace irred::cosmology {
// Positive differential opacity K=a*n_e*sigma_T in conformal Mpc^-1.
// Original binary64 knots define a linear law; no extrapolation or c factor.
struct SuppliedConformalOpacity {
  std::vector<double> eta_mpc, differential_opacity_mpc_inverse;
  std::string origin, law_id, exact_member_digest;
};
struct FiniteOpacityRequest {
  ThermalPhysicalModel model{}; // no other radiation or explicit species
  double initial_scale_factor = 0, observer_scale_factor = 0;
  SuppliedConformalOpacity opacity;
  std::vector<double> k_mpc_inverse;
  std::string source_origin;
};
struct FiniteOpacityPolicy {
  double state_absolute_tolerance = 1e-8, relative_tolerance = 3e-5;
  double channel_absolute_tolerance_mpc_inverse = 1e-11;
  double maximum_eta_step_mpc = 4, maximum_log_a_step = .03;
  std::size_t maximum_wavenumbers = 2, maximum_original_times = 1024;
  std::size_t maximum_fine_times = 4096, maximum_hierarchy = 192;
  std::size_t maximum_attempted_steps = 100000;
  std::size_t maximum_background_clock_calls = 1000000;
  // One attempt to solve the WHOLE two-stage perturbation system. Individual
  // 2x2 inversions/core factors have separate actual counters and write costs.
  std::size_t maximum_coupled_stage_solves = 200000;
  std::size_t maximum_destination_writes = 2000000000;
  std::size_t maximum_native_bytes = 64 * 1024 * 1024;
};
struct FiniteOpacityWork {
  std::size_t attempted_steps = 0, background_clock_calls = 0;
  std::size_t coupled_stage_solves = 0, tail_block_inversions = 0;
  std::size_t core_factorizations = 0, destination_writes = 0;
  // destination_writes is a charged logical buffer-update ledger, including
  // initialization/copy allowances. It is not a count of compiler/CPU stores;
  // scalar expression temporaries and shared thermal leaves have separate calls.
  // Space held for mandatory full-tail refusal receipts, not executed writes.
  std::size_t reserved_refusal_writes = 0;
  std::size_t refused_write_request = 0;
};
// Every native core coordinate is dimensionless. Velocities are theta/k.
// Order: delta_gamma, theta_gamma/k, sigma_gamma, G0,G1,G2,
//        delta_b,theta_b/k,delta_c,theta_c/k,phi.
using FiniteOpacityCore = std::array<long double, 11>;
struct FiniteOpacitySeed {
  std::array<double, 11> emitted_core{};
  std::array<long double, 11> absolute_cast_loss{};
  long double hamiltonian_residual = 0, zeta_residual = 0;
};
struct FiniteOpacityClockReceipt {
  long double initial_age_mpc = 0, observer_age_mpc = 0;
  long double initial_quadrature_difference_mpc = 0;
  long double observer_quadrature_difference_mpc = 0;
  long double initial_eta_mismatch_mpc = 0, observer_eta_mismatch_mpc = 0;
  // Differences are witnesses, not complete clock/source propagation errors.
  bool endpoint_clock_consistent = false;
};
struct FiniteOpacityBoundary {
  long double eta_i_mpc = 0, tau_i = 0, survival_i = 0;
  std::vector<long double> photon_monopole, photon_dipole_theta_over_k;
  std::vector<long double> omitted_temperature_absolute_bound;
  // Every F_l>=2 and every G_l is exactly zero in this compiled seed. Thus
  // B_l=w_i*(monopole*j_l+dipole*j_l') and initial E is exactly zero.
};
struct FiniteOpacityIdentity {
  std::string_view model_id, mode_id, method_id, arithmetic_id;
  FiniteOpacityRequest original;
  FiniteOpacityPolicy preparation_policy;
  std::array<ThermalScalarMapWitness, 4> mapping{};
  ThermalBackground background;
  long double lambda_retained = 0, lambda_getter_signed_loss = 0;
  long double critical_density_ev4 = 0, species_normalization_error = 0;
  double lambda_emitted = 0;
  FiniteOpacityClockReceipt clock;
  std::vector<FiniteOpacitySeed> seeds;
  FiniteOpacityBoundary boundary;
  FiniteOpacityWork preparation_work;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
};
struct FiniteOpacityNode {
  std::size_t eta_index = 0, k_index = 0;
  long double eta_mpc = 0, scale_factor = 0;
  FiniteOpacityCore core{};
  std::array<long double, 4> raw_channels{}, channel_cast_loss{};
  long double psi = 0, phi_prime_mpc_inverse = 0;
  long double hamiltonian_residual = 0, hamiltonian_term_scale = 0;
  long double normalized_hamiltonian_residual = 0;
  long double closure_defect = 0, actual_p_minus_shadow = 0;
  long double conditional_hcal_estimate = 0;
  long double survival = 0, visibility_mpc_inverse = 0;
};
struct FiniteOpacityAttemptReceipt {
  unsigned hierarchy = 0, time_refinement = 0;
  numerics::Status status = numerics::Status::invalid_input;
  std::vector<FiniteOpacityNode> nodes; // k-major, all explicit indices
  // Full finest reached hierarchy per k, including a refused final prefix.
  std::vector<long double> final_temperature_tail, final_polarization_tail;
  std::vector<FiniteOpacityCore> final_core;
  std::vector<long double> final_eta_mpc, final_scale_factor;
  std::size_t attempted_steps = 0, completed_steps = 0;
  std::size_t reached_wavenumbers = 0;
  long double maximum_stage_residual = 0, minimum_scaled_pivot = 1;
  long double maximum_clock_stage_residual = 0;
  long double maximum_observer_scale_factor_mismatch = 0;
  bool observer_scale_factor_consistent = false;
};
struct FiniteOpacityChannelDiagnostic {
  std::size_t eta_index = 0, k_index = 0;
  unsigned channel = 0;
  long double allocation = 0, time_difference = 0, hierarchy_difference = 0;
  long double previous_time_difference = 0, previous_hierarchy_difference = 0;
  long double measured_cast_loss = 0;
  // Missing integrated error owners stay absent, even after endpoint agreement.
  std::optional<long double> common_background_clock_error;
  std::optional<long double> arithmetic_linear_error, source_grid_error;
};
struct FiniteOpacitySourceResult {
  FiniteOpacitySourceResult() = default;
  FiniteOpacitySourceResult(const FiniteOpacitySourceResult &) = delete;
  FiniteOpacitySourceResult &operator=(const FiniteOpacitySourceResult &) = delete;
  FiniteOpacitySourceResult(FiniteOpacitySourceResult &&) noexcept = default;
  FiniteOpacitySourceResult &operator=(FiniteOpacitySourceResult &&) noexcept = default;
  numerics::Status status = numerics::Status::invalid_input;
  FiniteOpacityPolicy policy;
  std::shared_ptr<const FiniteOpacityIdentity> identity;
  // Borrowed from the SAME immutable identity, preserved on every evaluation
  // refusal and after producer destruction. No second vector owner is made.
  const FiniteOpacityBoundary *boundary() const noexcept {
    return identity ? &identity->boundary : nullptr;
  }
  // Exactly one COMPLETE computed source, never a refused prefix. It is raw
  // engineering output while admission is false; projection cannot infer this
  // producer's physical/error admission from the strings on its source.
  std::optional<projection::ContinuousCmbSource> source;
  bool source_numerically_admitted = false;
  std::vector<FiniteOpacityAttemptReceipt> attempts;
  std::vector<FiniteOpacityChannelDiagnostic> diagnostics;
  FiniteOpacityWork work;
  std::optional<std::size_t> peak_owned_payload_bound;
};
class FiniteOpacitySourceProducer {
public:
  FiniteOpacitySourceProducer() = default;
  FiniteOpacitySourceProducer(const FiniteOpacitySourceProducer &) = delete;
  FiniteOpacitySourceProducer &operator=(const FiniteOpacitySourceProducer &) = delete;
  FiniteOpacitySourceProducer(FiniteOpacitySourceProducer &&) noexcept;
  FiniteOpacitySourceProducer &operator=(FiniteOpacitySourceProducer &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const FiniteOpacityIdentity *identity() const noexcept { return identity_.get(); }
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  FiniteOpacitySourceResult produce(FiniteOpacityPolicy = {}) const;
private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::shared_ptr<const FiniteOpacityIdentity> identity_;
  friend FiniteOpacitySourceProducer prepare_finite_opacity_source(
      FiniteOpacityRequest &&, FiniteOpacityPolicy);
};
FiniteOpacitySourceProducer prepare_finite_opacity_source(
    FiniteOpacityRequest &&, FiniteOpacityPolicy = {});
inline constexpr std::string_view finite_opacity_model_id =
    "GR/flat-photon-cold-baryon-CDM-smooth-Lambda-supplied-finite-opacity-finite-adiabatic-boundary/v1";
inline constexpr std::string_view finite_opacity_mode_id =
    "conformal-Newtonian-finite-start-entropy-free-local-zeta-one-zero-higher-photon-moments/v1";
inline constexpr std::string_view finite_opacity_method_id =
    "full-temperature-polarization-hierarchy-RadauIIA2-block-Schur-time-L-refinement-raw-split-source/v1";
inline constexpr std::string_view finite_opacity_arithmetic_id =
    "emitted-binary64-input-strict-wide-nearest-shared-thermal-conformal-MB-moments/v1";
} // namespace irred::cosmology
