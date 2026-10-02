#pragma once
#include "irred/hydrogen_equilibrium.hpp"
#include "irred/thermal_neutrino.hpp"
namespace irred::cosmology {
struct PureHydrogenRequest {
  ThermalPhysicalModel model;
  double initial_redshift = 1600, late_redshift = 300;
};
struct PureHydrogenPolicy {
  double absolute_x_tolerance = 1e-8, relative_x_tolerance = 1e-6;
  double absolute_tau_tolerance = 2e-7, relative_tau_tolerance = 1e-6;
  double absolute_root_tolerance = 2e-3;
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
                          hydrogen_drag_depth = 2;
struct HydrogenHistoryRow {
  double redshift = 0;
  HydrogenHistoryValue electron_fraction, drag_depth;
};
struct HydrogenHistoryBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested = 0;
  std::vector<HydrogenHistoryRow> rows;
};
struct HydrogenHistoryWork {
  std::size_t background_evaluations = 0, momentum_callbacks = 0,
              rhs_evaluations = 0, equilibrium_solves = 0;
  std::size_t total() const noexcept {
    return background_evaluations + momentum_callbacks + rhs_evaluations +
           equilibrium_solves;
  }
};
// Compiled pure-H F=1 Peebles3level, prescribed Tm=Tr, not evolved temperature.
// The optical depth and unit-depth root use the explicit late endpoint >=300;
// they are not a full cosmological drag epoch or a replacement for
// supplieddrag.
class PureHydrogenHistory {
public:
  PureHydrogenHistory() = default;
  PureHydrogenHistory(const PureHydrogenHistory &) = default;
  PureHydrogenHistory &operator=(const PureHydrogenHistory &) = default;
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
  pure_hydrogen_payload_bound(std::size_t, std::size_t, std::size_t) noexcept;
};
PureHydrogenHistory prepare_pure_hydrogen_history(const PureHydrogenRequest &,
                                                  PureHydrogenPolicy = {});
std::optional<std::size_t>
pure_hydrogen_payload_bound(std::size_t fine_intervals,
                            std::size_t output_points,
                            std::size_t species) noexcept;
inline constexpr std::string_view pure_hydrogen_history_id =
    "pure-H-Peebles3level-F1-prescribed-Tm-equals-Tr-truncated-drag/v1";
} // namespace irred::cosmology
