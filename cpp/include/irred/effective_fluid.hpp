#pragma once
#include "irred/thermal_observables.hpp"
namespace irred::cosmology {
// Independently declared nominal fluid; not a source-normalized scalar EDE.
// reference contains the complete once-emitted thermal mapping, unchanged.
struct EffectiveFluidRequest {
  ThermalFlatModel reference;
  double density_at_transition = 0, transition_scale = 0, ruler_endpoint = 0;
  std::string source_origin, endpoint_origin;
};
struct EffectiveFluidState {
  numerics::Status status = numerics::Status::invalid_input;
  long double density = 0, equation_of_state = 0, log_scale_derivative = 0;
  long double density_estimate = 0, equation_of_state_estimate = 0;
  long double derivative_estimate = 0;
};
inline constexpr unsigned effective_fluid_ruler = 4;
struct EffectiveFluidBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested_outputs = 0;
  std::optional<EarlyLateValue> ruler;
  std::vector<ThermalBackgroundRow> rows;
  std::size_t callbacks = 0, outer_callbacks = 0, momentum_callbacks = 0;
};
class EffectiveFluidBackground {
public:
  EffectiveFluidBackground() = default;
  EffectiveFluidBackground(const EffectiveFluidBackground &) = default;
  EffectiveFluidBackground &operator=(const EffectiveFluidBackground &);
  EffectiveFluidBackground(EffectiveFluidBackground &&) noexcept;
  EffectiveFluidBackground &operator=(EffectiveFluidBackground &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const EffectiveFluidRequest &request() const noexcept { return request_; }
  const ThermalFlatModel &source() const noexcept { return request_.reference; }
  const ThermalBackground &reference_background() const noexcept { return reference_; }
  std::size_t preparation_callbacks() const noexcept { return reference_.preparation_callbacks(); }
  std::optional<long double> omega_lambda() const noexcept;
  EffectiveFluidState fluid_state(long double a) const noexcept;
  ThermalScaledExpansion scaled_expansion(long double a, ThermalPolicy = {}) const;
  EffectiveFluidBatch evaluate(std::span<const double> scale_factors,
                              unsigned requested_outputs,
                              ThermalObservablePolicy = {}) const;
private:
  numerics::Status status_ = numerics::Status::invalid_input;
  EffectiveFluidRequest request_{};
  ThermalBackground reference_;
  long double today_density_ = 0, today_estimate_ = 0, transition_power_ = 0;
  long double lambda_ = 0;
  friend EffectiveFluidBackground prepare_effective_fluid(
      const EffectiveFluidRequest &, ThermalObservablePolicy);
};
EffectiveFluidBackground prepare_effective_fluid(const EffectiveFluidRequest &,
                                                ThermalObservablePolicy = {});
// PREP/EVAL/move-return scope: this retained owner and operation temporaries/
// rows. Excludes borrowed inputs, extra caller-owned copies, arbitrary copy
// assignment's old-target/replacement coexistence, stack, allocator/RSS.
// Copy operations have no policy parameter; callers own their whole lifetime.
std::optional<std::size_t> effective_fluid_payload_bound(
    std::size_t points, std::size_t species, std::size_t origin_bytes) noexcept;
inline constexpr std::string_view effective_fluid_model_id =
    "independent-compensated-density-nine-halves-thermal-FD-background/v1";
inline constexpr std::string_view effective_fluid_ruler_id =
    "independent-effective-fluid-matched-finite-photon-baryon-ruler/v1";
} // namespace irred::cosmology
