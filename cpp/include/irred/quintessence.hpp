#pragma once
#include "irred/numerics.hpp"
#include <string_view>
namespace irred::cosmology {
// Anchor a=1; reduced Planck mass; positive exponential potential or its
// exactly zero boundary. Dust and radiation are separately conserved in GR.
struct ExponentialQuintessence {
  double lambda = 0;
  double h_anchor_km_s_mpc = 70;
  double signed_kinetic_fraction_root = 0;
  double potential_fraction = .7;
  double radiation_fraction = 0;
};
// Tolerances allocate local step diagnostics and the Friedmann residual.
// They do not promise final observable accuracy after dynamical amplification.
struct QuintessencePolicy {
  double absolute_tolerance = 1e-10, relative_tolerance = 1e-9;
  std::size_t maximum_callbacks_per_point = 100000,
              maximum_total_callbacks = 1000000;
  std::size_t maximum_points = 4096, maximum_native_bytes = 16 * 1024 * 1024;
  unsigned maximum_halvings = 30;
};
struct QuintessenceValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  // Local step-doubling/arithmetic sums transformed at the endpoint; no
  // propagated global state error or certified observable bound is supplied.
  double absolute_error_estimate = 0;
};
struct QuintessenceRow {
  double scale_factor = 0;
  numerics::Status status = numerics::Status::invalid_input;
  QuintessenceValue e, h_km_s_mpc, w_phi, omega_phi, omega_m, omega_r;
  QuintessenceValue dm_mpc, dl_mpc;
  double constraint_residual = 0;
  std::size_t callbacks = 0, accepted_steps = 0, rejected_steps = 0;
};
struct QuintessenceBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::vector<QuintessenceRow> rows;
  std::size_t callbacks = 0;
};
class QuintessenceBackground {
public:
  QuintessenceBackground() = default;
  numerics::Status status() const noexcept { return status_; }
  const ExponentialQuintessence &specification() const noexcept { return source_; }
  QuintessenceBatch evaluate(std::span<const double> scale_factors,
                            bool distances = false,
                            QuintessencePolicy = {}) const;
private:
  ExponentialQuintessence source_{};
  numerics::Status status_ = numerics::Status::invalid_input;
  friend QuintessenceBackground prepare_quintessence(ExponentialQuintessence);
};
QuintessenceBackground prepare_quintessence(ExponentialQuintessence);
std::optional<std::size_t> quintessence_payload_bound(std::size_t) noexcept;
inline constexpr std::string_view quintessence_model_id =
    "GR/flat-canonical-exponential-quintessence-dust-radiation-IVP/v1";
} // namespace irred::cosmology
