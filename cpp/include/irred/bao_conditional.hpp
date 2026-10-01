#pragma once
#include "irred/bao.hpp"
#include "irred/early_late.hpp"
namespace irred::bao {
struct ConditionalDensityPolicy {
  cosmology::EarlyLatePolicy predictions;
  size_t maximum_models = 0, maximum_queries = 0, maximum_string_bytes = 0,
         maximum_native_bytes = 0, maximum_total_callbacks = 0;
  double maximum_forward_sensitivity = 0;
  double maximum_projection_log_density_error = 1e-8;
  numerics::Arithmetic arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  std::uint32_t requested = 0; // existing Output bits
};
struct ConditionalDensitySlot {
  cosmology::SoundHorizonRequest source;
  OutputState predictions_state, residuals_state, density_state;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::optional<statistics::GaussianResult> result;
  std::vector<double> predictions, residuals;
  std::optional<double> projection_log_density_error_estimate;
  size_t callbacks = 0;
};
struct ConditionalDensityBatch {
  statistics::DensityStatus status = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<ConditionalDensitySlot> slots;
  size_t callbacks = 0;
};
// Simultaneous returned payload and serial scratch, excluding retained source,
// borrowed inputs, allocator bookkeeping, stack frames and RSS.
// origin_bytes is the sum of max(32, origin.size()+1) for copied requests.
std::optional<size_t>
conditional_density_payload_bound(size_t models, size_t rows,
                                  size_t origin_bytes,
                                  unsigned requested) noexcept;
inline constexpr std::string_view conditional_density_id =
    "BAO/flat-early-late-supplied-drag-normalized-ratio-density/v1";
} // namespace irred::bao
