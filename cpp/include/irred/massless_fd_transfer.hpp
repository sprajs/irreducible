#pragma once
#include "irred/thermal_neutrino.hpp"
#include <array>

namespace irred::cosmology {
inline constexpr unsigned massless_fd_comoving_cdm = 1,
                          massless_fd_spatial_potential = 2,
                          massless_fd_lapse_potential = 4;
// Fixed finite-domain trial admission. Smaller tolerances/caps may be
// requested; the accepted source contract's tolerances/caps cannot be widened
// silently.
struct MasslessFDTransferPolicy {
  double absolute_tolerance = 1e-7, relative_tolerance = 1e-4;
  double maximum_log_step = .04, maximum_phase_step = .04;
  double maximum_constraint_residual = 1e-7;
  std::size_t maximum_points = 16;
  std::size_t maximum_rhs_per_point = 1000000, maximum_total_rhs = 4000000;
  std::size_t maximum_scalar_updates_per_point = 128000000,
              maximum_total_scalar_updates = 512000000;
  std::size_t maximum_background_queries = 4000000;
  std::size_t maximum_age_quadrature_evaluations = 200000;
  std::size_t maximum_native_bytes = 16 * 1024 * 1024;
  ThermalPolicy thermal;
};
struct MasslessFDTransferWork {
  std::size_t rhs = 0, scalar_updates = 0, background_queries = 0;
  std::size_t age_quadrature_evaluations = 0, momentum_callbacks = 0;
};
struct MasslessFDTransferValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double absolute_error_estimate = 0;
  double time_refinement = 0, initial_refinement = 0, hierarchy_refinement = 0,
         background_age_sensitivity = 0, arithmetic_cast_sensitivity = 0;
};
// All continuous leading/projected states survive arithmetic admission. The
// signed Newtonian delta_c witness is distinct from tiny comoving Delta_c.
struct MasslessFDInitialWitness {
  numerics::Status status = numerics::Status::invalid_input;
  double scale_factor = 0;
  long double eta_mpc = 0, eta_error_mpc = 0, hcal_mpc_inverse = 0;
  long double leading_phi = 0, leading_psi = 0, leading_delta_c = 0,
              leading_delta_r = 0, leading_theta = 0, leading_sigma = 0;
  long double projected_phi = 0, projected_psi = 0, projected_vc = 0,
              projected_delta = 0, projected_phi_n = 0, scaled_shear = 0;
  long double lapse_correction = 0, velocity_correction = 0;
  long double hamiltonian_before = 0, hamiltonian_after = 0;
  long double radiation_fraction = 0, closure_defect = 0;
};
struct MasslessFDRunWitness {
  numerics::Status status = numerics::Status::invalid_input;
  bool endpoint_available = false;
  bool lapse_available = false;
  unsigned hierarchy_l = 0, time_divisor = 0, control_kind = 0;
  int control_sign = 0;
  MasslessFDInitialWitness initial;
  // Delta and phi survive a failed final background query; psi is present only
  // when lapse_available binds the metric epoch to this state epoch.
  std::array<long double, 3> endpoint{};
  long double endpoint_scale_factor = 0, metric_epoch_scale_factor = 0;
  long double endpoint_scaled_shear = 0;
  long double endpoint_eta_mpc = 0;
  long double maximum_stage_phase_bound = 0, maximum_phase_increment_upper = 0;
  double maximum_constraint_residual = 0;
  MasslessFDTransferWork work;
};
struct MasslessFDTransferRow {
  double wavenumber_mpc_inverse = 0;
  MasslessFDTransferValue comoving_cdm, spatial_potential, lapse_potential;
  double maximum_constraint_residual = 0;
  double maximum_closure_defect = 0, maximum_background_identity_defect = 0;
  long double conformal_age_mpc = 0, conformal_age_error_mpc = 0;
  std::array<MasslessFDInitialWitness, 3> initial_states{};
  // Eight anchors, eight matched numerical-input perturbations, one arithmetic
  // control. Failed trials retain their initial/endpoint/status/counter record.
  std::array<MasslessFDRunWitness, 17> attempts{};
  unsigned attempts_recorded = 0;
  MasslessFDTransferWork work;
};
struct MasslessFDTransferBatch {
  numerics::Status status = numerics::Status::invalid_input;
  double scale_factor = 0;
  unsigned requested_outputs = 0;
  std::vector<MasslessFDTransferRow> rows;
  MasslessFDTransferWork work;
};
class MasslessFDTransfer {
public:
  MasslessFDTransfer() = default;
  MasslessFDTransfer(const MasslessFDTransfer &) = default;
  MasslessFDTransfer &operator=(const MasslessFDTransfer &) = default;
  MasslessFDTransfer(MasslessFDTransfer &&) noexcept;
  MasslessFDTransfer &operator=(MasslessFDTransfer &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const ThermalBackground *background() const noexcept {
    return background_ ? &*background_ : nullptr;
  }
  double initial_scale_factor() const noexcept { return initial_; }
  // Preserve the caller's k order and repeated modes. Each field is a signed
  // transfer per unit asymptotic zeta, not a spectrum or a supplied amplitude.
  MasslessFDTransferBatch evaluate(std::span<const double> ordered_k,
                                   double common_a, unsigned requested_outputs,
                                   MasslessFDTransferPolicy = {}) const;

private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<ThermalBackground> background_;
  double initial_ = 0;
  friend MasslessFDTransfer
  prepare_massless_fd_transfer(const ThermalBackground &, double);
};
MasslessFDTransfer prepare_massless_fd_transfer(const ThermalBackground &,
                                                double initial_a);
std::optional<std::size_t>
massless_fd_transfer_payload_bound(std::size_t points,
                                   std::size_t species) noexcept;
inline constexpr std::string_view massless_fd_transfer_model_id =
    "GR/flat-pure-explicit-massless-FD-CDM-lambda-unit-zeta/v1";
inline constexpr std::string_view massless_fd_transfer_method_id =
    "projected-leading-adiabatic-RK4-log-a-L32-64-128-start-time-refinement/v1";
inline constexpr std::string_view massless_fd_transfer_arithmetic_id =
    "strict-wide-scaled-shear-comoving-state-unit-zeta/v1";
} // namespace irred::cosmology
