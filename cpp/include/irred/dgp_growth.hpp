#pragma once
#include "irred/numerics.hpp"
#include <string_view>
namespace irred::cosmology {
inline constexpr std::string_view dgp_growth_id =
    "DGP/flat-self-accelerating-dust-quasistatic-growing-mode/v1";
inline constexpr unsigned dgp_e = 1, dgp_h = 2, dgp_omega_m = 4, dgp_mu = 8,
                          dgp_d = 16, dgp_f = 32;
struct DGPParameters {
  double omega_m0 = .3, h0_km_s_mpc = 70;
};
struct DGPPolicy {
  double relative_tolerance = 1e-8;
  std::size_t maximum_points = 4096, maximum_native_bytes = 16 * 1024 * 1024;
  std::size_t maximum_steps_per_point = 20000, maximum_total_steps = 200000;
};
struct DGPValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double absolute_error_estimate = 0;
};
struct DGPRow {
  double scale_factor = 0;
  DGPValue e, h, omega_m, mu, d, f;
  std::size_t trials = 0, rejected_trials = 0, callbacks = 0;
};
struct DGPBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested = 0;
  std::vector<DGPRow> rows;
  std::size_t trials = 0, rejected_trials = 0, callbacks = 0;
};
class DGPGrowth {
public:
  DGPGrowth() = default;
  DGPGrowth(const DGPGrowth &) = default;
  DGPGrowth &operator=(const DGPGrowth &) = default;
  DGPGrowth(DGPGrowth &&) noexcept;
  DGPGrowth &operator=(DGPGrowth &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const DGPParameters *parameters() const noexcept {
    return parameters_ ? &*parameters_ : nullptr;
  }
  DGPBatch evaluate(std::span<const double> scale_factors, unsigned outputs,
                    DGPPolicy = {}) const;

private:
  std::optional<DGPParameters> parameters_;
  numerics::Status status_ = numerics::Status::invalid_input;
  friend DGPGrowth prepare_dgp_growth(DGPParameters);
};
DGPGrowth prepare_dgp_growth(DGPParameters);
std::optional<std::size_t> dgp_payload_bound(std::size_t points) noexcept;
} // namespace irred::cosmology
