#pragma once
#include "irred/early_late.hpp"
#include "irred/thermal_neutrino.hpp"
namespace irred::cosmology {
struct ThermalObservableRequest {
  ThermalPhysicalModel model;
  double z_drag;
  std::string drag_origin, source_origin;
};
// The nine EarlyLateOutput coordinates retain their units and meanings.
// This independent bit requests the batch's single retained ruler directly.
inline constexpr unsigned thermal_ruler_mask = 1u << early_late_output_count;
struct ThermalObservablePolicy {
  double absolute_tolerance_mpc = 1e-8, relative_tolerance = 2e-10;
  double absolute_tolerance_ratio = 1e-10, relative_tolerance_ratio = 5e-10;
  std::size_t maximum_callbacks_per_point = 20000000;
  std::size_t maximum_total_callbacks = 100000000;
  unsigned maximum_depth = 30;
  std::size_t maximum_points = 4096, maximum_native_bytes = 16 * 1024 * 1024;
  ThermalPolicy thermal = [] {
    auto p = ThermalPolicy{};
    // A complete nested distance/ruler batch needs more aggregate work than
    // one background-only batch; its separate moment ceiling stays explicit.
    p.maximum_total_callbacks = 100000000;
    return p;
  }();
};
struct ThermalObservableBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested_outputs = 0;
  std::optional<EarlyLateValue> ruler;
  std::vector<EarlyLateRow> rows;
  std::size_t callbacks = 0, outer_callbacks = 0, momentum_callbacks = 0;
};
class ThermalObservables {
public:
  ThermalObservables() = default;
  ThermalObservables(const ThermalObservables &) = default;
  ThermalObservables &operator=(const ThermalObservables &) = default;
  ThermalObservables(ThermalObservables &&) noexcept;
  ThermalObservables &operator=(ThermalObservables &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const ThermalObservableRequest &source() const noexcept { return source_; }
  const ThermalBackground &background() const noexcept { return background_; }
  ThermalObservableBatch evaluate(std::span<const double> redshifts,
                                  unsigned requested_outputs,
                                  ThermalObservablePolicy = {}) const;

private:
  numerics::Status status_ = numerics::Status::invalid_input;
  ThermalObservableRequest source_{};
  ThermalBackground background_;
  friend ThermalObservables
  prepare_thermal_observables(const ThermalObservableRequest &,
                              ThermalObservablePolicy);
};
ThermalObservables prepare_thermal_observables(const ThermalObservableRequest &,
                                               ThermalObservablePolicy = {});
// Includes retained physical/mapped owners, mapping/prepare temporaries, row
// capacities, origin copies and nested moment scratch. Excludes allocator
// metadata, recursion stack and RSS; borrowed input owners are excluded.
std::optional<std::size_t>
thermal_observables_payload_bound(std::size_t points, std::size_t species,
                                  std::size_t origin_bytes) noexcept;
inline constexpr std::string_view thermal_observables_equation_id =
    "flat-thermal-FD-log-redshift-distance-supplied-drag-ruler/v1";
inline constexpr std::string_view thermal_physical_mapping_id =
    "physical-omega-h2-explicit-Kelvin-temperature-state-weight-SI/v1";
} // namespace irred::cosmology
