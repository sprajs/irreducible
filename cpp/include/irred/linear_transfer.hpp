#pragma once
#include "irred/thermal_neutrino.hpp"
#include <array>
namespace irred::cosmology {
inline constexpr unsigned transfer_comoving_cdm = 1, transfer_metric = 2;
struct TransferPolicy {
  double absolute_tolerance = 1e-8, relative_tolerance = 3e-5;
  double maximum_log_step = .08, maximum_constraint_residual = 1e-6;
  std::size_t maximum_points = 4096, maximum_total_callbacks = 20000000;
  std::size_t maximum_callbacks_per_point = 4000000;
  std::size_t maximum_native_bytes = 16 * 1024 * 1024;
};
struct TransferValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double absolute_error_estimate = 0, time_refinement = 0,
         initial_refinement = 0;
};
struct TransferRow {
  double wavenumber_mpc_inverse = 0;
  TransferValue comoving_cdm, metric;
  double maximum_constraint_residual = 0;
  std::size_t callbacks = 0;
};
struct TransferBatch {
  numerics::Status status = numerics::Status::invalid_input;
  double scale_factor = 0;
  unsigned requested_outputs = 0;
  std::vector<TransferRow> rows;
  std::size_t callbacks = 0;
};
struct PrimordialBand {
  double amplitude, spectral_index, pivot_mpc_inverse;
  double minimum_mpc_inverse, maximum_mpc_inverse;
};
struct BandVariancePolicy {
  TransferPolicy transfer;
  double absolute_tolerance = 1e-10, relative_tolerance = 1e-3;
  unsigned base_panels = 8; // n/2n/4n composite Simpson; n must be even.
};
struct BandVariance {
  numerics::Status status = numerics::Status::invalid_input;
  PrimordialBand primordial{};
  double scale_factor = 0, radius_mpc = 0;
  std::optional<double> variance, sigma;
  // The first diagnostic has variance units; sigma has its own propagation.
  double sigma_absolute_error_estimate = 0;
  double absolute_error_estimate = 0, window_refinement = 0,
         transfer_error_estimate = 0;
  std::size_t callbacks = 0;
  // Exactly zero outside the supplied primordial band; no ordinary sigma8 tail.
  bool primordial_support_is_finite_band = true;
};
class PerfectFluidTransfer {
public:
  PerfectFluidTransfer() = default;
  PerfectFluidTransfer(const PerfectFluidTransfer &) = default;
  PerfectFluidTransfer &operator=(const PerfectFluidTransfer &) = default;
  PerfectFluidTransfer(PerfectFluidTransfer &&) noexcept;
  PerfectFluidTransfer &operator=(PerfectFluidTransfer &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const ThermalBackground *background() const noexcept {
    return background_ ? &*background_ : nullptr;
  }
  double initial_scale_factor() const noexcept { return initial_; }
  // Signed unit-zeta transfer; caller k order is retained, no sorting.
  TransferBatch evaluate(std::span<const double> k_mpc_inverse,
                         double scale_factor, unsigned outputs,
                         TransferPolicy = {}) const;
  BandVariance band_variance(const PrimordialBand &, double scale_factor,
                             double radius_mpc, BandVariancePolicy = {}) const;
  BandVariance sigma8_band(const PrimordialBand &, double scale_factor,
                           BandVariancePolicy = {}) const;

private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<ThermalBackground> background_;
  double initial_ = 0;
  friend PerfectFluidTransfer
  prepare_perfect_fluid_transfer(const ThermalBackground &, double);
};
PerfectFluidTransfer
prepare_perfect_fluid_transfer(const ThermalBackground &,
                               double initial_scale_factor);
std::optional<std::size_t> transfer_payload_bound(std::size_t points) noexcept;
inline constexpr std::string_view perfect_fluid_transfer_model_id =
    "GR/flat-self-interacting-perfect-radiation-CDM-lambda-unit-zeta/v1";
inline constexpr std::string_view perfect_fluid_transfer_method_id =
    "constraint-consistent-adiabatic-series-phase-resolved-RK4-start-time-"
    "refinement/v1";
inline constexpr std::string_view perfect_fluid_transfer_arithmetic_id =
    "strict-wide-log-a-cancellation-safe-comoving-state/v1";
inline constexpr std::string_view primordial_band_variance_id =
    "unit-zeta-power-law-compact-k-support-CDM-top-hat/v1";
} // namespace irred::cosmology
