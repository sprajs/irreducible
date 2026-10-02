#pragma once
#include "irred/bao.hpp"
#include "irred/thermal_observables.hpp"
namespace irred::bao {
struct ThermalDensityPolicy {
  cosmology::ThermalObservablePolicy predictions;
  size_t maximum_models = 0, maximum_queries = 0, maximum_string_bytes = 0,
         maximum_native_bytes = 0, maximum_total_callbacks = 0;
  double maximum_forward_sensitivity = 0;
  double maximum_projection_log_density_error = 1e-8;
  numerics::Arithmetic arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  std::uint32_t requested = 0; // existing Output bits
};
struct ThermalDensitySlot {
  cosmology::ThermalObservableRequest source;
  OutputState predictions_state, residuals_state, density_state;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  numerics::Status preparation_status = numerics::Status::invalid_input;
  std::optional<statistics::GaussianResult> result;
  std::vector<double> predictions, residuals;
  std::optional<double> projection_log_density_error_estimate;
  size_t callbacks = 0, preparation_callbacks = 0, outer_callbacks = 0,
         momentum_callbacks = 0;
};
struct ThermalDensityBatch {
  statistics::DensityStatus status = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<ThermalDensitySlot> slots;
  size_t callbacks = 0, preparation_callbacks = 0, outer_callbacks = 0,
         momentum_callbacks = 0;
};
// Simultaneous returned payload and serial scratch, excluding retained source,
// borrowed inputs, allocator bookkeeping, stack frames and RSS.
// origin_bytes sums max(32, size()+1) for both origins per request.
// species counts all copied slot species; maximum_species bounds serial
// provider.
std::optional<size_t>
thermal_density_payload_bound(size_t models, size_t rows, size_t origin_bytes,
                              size_t species, size_t maximum_species,
                              unsigned requested) noexcept;
inline constexpr std::string_view thermal_density_id =
    "BAO/flat-thermal-FD-supplied-drag-normalized-ratio-density/v1";
} // namespace irred::bao
