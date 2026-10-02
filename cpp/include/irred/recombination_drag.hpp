#pragma once
#include "irred/hydrogen_equilibrium.hpp"
#include "irred/thermal_neutrino.hpp"
namespace irred::cosmology {
enum class HydrogenTemperatureModel {
  prescribed_radiation,
  evolved_compton_adiabatic
};
struct PureHydrogenRequest {
  ThermalPhysicalModel model;
  double initial_redshift = 1600, late_redshift = 300;
  HydrogenTemperatureModel temperature_model =
      HydrogenTemperatureModel::prescribed_radiation;
};
struct PureHydrogenPolicy {
  double absolute_x_tolerance = 1e-8, relative_x_tolerance = 1e-6;
  double absolute_tau_tolerance = 2e-7, relative_tau_tolerance = 1e-6;
  double absolute_root_tolerance = 2e-3;
  double absolute_temperature_tolerance_kelvin = 2e-5,
         relative_temperature_tolerance = 2e-7;
  double absolute_thomson_tau_tolerance = 2e-7,
         relative_thomson_tau_tolerance = 1e-6;
  double absolute_opacity_tolerance = 2e-9, relative_opacity_tolerance = 2e-6;
  double absolute_visibility_tolerance = 2e-9,
         relative_visibility_tolerance = 3e-6;
  double absolute_survival_tolerance = 0, relative_survival_tolerance = 1e-5;
  std::size_t base_intervals = 8192, maximum_fine_intervals = 65536,
              maximum_total_work = 4000000, maximum_output_points = 4096,
              maximum_native_bytes = 16 * 1024 * 1024;
  ThermalPolicy thermal;
};
struct HydrogenHistoryValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double absolute_error_estimate = 0;
};
inline constexpr unsigned hydrogen_electron_fraction = 1,
                          hydrogen_drag_depth = 2,
                          hydrogen_matter_temperature = 4,
                          hydrogen_thomson_depth = 8,
                          hydrogen_thomson_opacity = 16,
                          hydrogen_visibility = 32, hydrogen_survival = 64;
struct HydrogenHistoryRow {
  double redshift = 0;
  HydrogenHistoryValue electron_fraction, drag_depth;
  HydrogenHistoryValue matter_temperature_kelvin, thomson_depth,
      thomson_opacity_per_redshift, visibility_per_redshift,
      finite_endpoint_survival;
};
struct HydrogenHistoryBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested = 0;
  std::vector<HydrogenHistoryRow> rows;
};
struct HydrogenHistoryWork {
  std::size_t background_evaluations = 0, momentum_callbacks = 0,
              rhs_evaluations = 0, equilibrium_solves = 0,
              temperature_rhs_evaluations = 0;
  std::size_t total() const noexcept {
    return background_evaluations + momentum_callbacks + rhs_evaluations +
           equilibrium_solves + temperature_rhs_evaluations;
  }
};
// Compiled pure-H F=1 Peebles3level; explicitly selected prescribed or evolved
// Tm. The optical depth and unit-depth root use the explicit late endpoint
// >=300; they are not a full cosmological drag epoch or a replacement for
// supplieddrag.
class PureHydrogenHistory {
public:
  PureHydrogenHistory() = default;
  PureHydrogenHistory(const PureHydrogenHistory &) = default;
  PureHydrogenHistory &operator=(const PureHydrogenHistory &);
  PureHydrogenHistory(PureHydrogenHistory &&) noexcept;
  PureHydrogenHistory &operator=(PureHydrogenHistory &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const PureHydrogenRequest *source() const noexcept {
    return source_ ? &*source_ : nullptr;
  }
  const ThermalBackground *background() const noexcept {
    return source_ ? &background_ : nullptr;
  }
  HydrogenHistoryWork work() const noexcept { return work_; }
  std::string_view model_identity() const noexcept;
  std::string_view method_identity() const noexcept;
  HydrogenHistoryValue conditional_unit_depth_redshift() const noexcept {
    return root_;
  }
  HydrogenHistoryBatch
  evaluate(std::span<const double> redshifts, unsigned outputs,
           std::size_t maximum_output_points = 4096,
           std::size_t maximum_native_bytes = 16 * 1024 * 1024) const;

private:
  struct Node {
    long double x = 0, x_error = 0, tau = 0, tau_error = 0;
  };
  struct ThermalNode {
    long double temperature = 0, temperature_error = 0, thomson_tau = 0,
                thomson_tau_error = 0, opacity = 0, opacity_error = 0;
  };
  std::vector<ThermalNode> thermal_nodes_;
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<PureHydrogenRequest> source_;
  ThermalBackground background_;
  PureHydrogenPolicy policy_;
  std::vector<Node> nodes_;
  HydrogenHistoryWork work_;
  HydrogenHistoryValue root_;
  friend PureHydrogenHistory
  prepare_pure_hydrogen_history(const PureHydrogenRequest &,
                                PureHydrogenPolicy);
  friend std::optional<std::size_t>
  pure_hydrogen_payload_bound(std::size_t, std::size_t, std::size_t,
                              HydrogenTemperatureModel) noexcept;
};
PureHydrogenHistory prepare_pure_hydrogen_history(const PureHydrogenRequest &,
                                                  PureHydrogenPolicy = {});
std::optional<std::size_t> pure_hydrogen_payload_bound(
    std::size_t fine_intervals, std::size_t output_points, std::size_t species,
    HydrogenTemperatureModel =
        HydrogenTemperatureModel::prescribed_radiation) noexcept;
inline constexpr std::string_view pure_hydrogen_history_id =
    "pure-H-Peebles3level-F1-prescribed-Tm-equals-Tr-truncated-drag/v1";
inline constexpr std::string_view evolved_hydrogen_history_id =
    "pure-H-Peebles3level-F1-SSS1999-Tm-rates-Compton-adiabatic-finite-"
    "endpoint-Thomson/v1";
} // namespace irred::cosmology
