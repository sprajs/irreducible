#pragma once
#include "irred/numerics.hpp"
#include <string>
#include <string_view>
namespace irred::cosmology {
inline constexpr std::string_view helium_escape_history_model_id =
    "HeI-HyRec2020-source-effective-escape-fixed-alpha-me-"
    "supplied-H-neutral-H-nuclei-TrTm/v1";
inline constexpr std::string_view helium_escape_history_method_id =
    "scalar-BE-ln-a-driver-knot-split-Richardson2-refinement3/v1";
inline constexpr std::string_view helium_escape_history_arithmetic_id =
    "longdouble64-cpu-nearest-strict-fp/v1";
inline constexpr std::string_view helium_escape_history_assets_id =
    "CLASS0ceb-HyRec2020-helium-cc113cdc-SI-fixed-native-sigmaT/v1";
enum class HeliumEscapeDriverRole {
  supplied_cosmological_history,
  synthetic_prescribed_bath
};
struct HeliumEscapeDriverKnot {
  double redshift = 0;
  double hubble_per_second = 0;
  double hydrogen_nuclei_per_cubic_metre = 0;
  double radiation_temperature_kelvin = 0;
  double hydrogen_neutral_fraction = 0;
};
// Borrowed source; prepare acquires exactly one owned driver. No abundance,
// background, H kinetics or temperature law is inferred from these columns.
struct HeliumEscapeHistoryRequest {
  std::span<const HeliumEscapeDriverKnot> knots;
  double helium_to_hydrogen_nuclei_ratio = 0;
  double initial_helium_singly_ionized_fraction = 0;
  std::string_view producer_identity;
  HeliumEscapeDriverRole role = HeliumEscapeDriverRole::supplied_cosmological_history;
};
struct HeliumEscapeDriver {
  std::vector<HeliumEscapeDriverKnot> knots;
  double helium_to_hydrogen_nuclei_ratio = 0;
  double initial_helium_singly_ionized_fraction = 0;
  std::string producer_identity;
  HeliumEscapeDriverRole role = HeliumEscapeDriverRole::supplied_cosmological_history;
};
struct HeliumEscapeHistoryPolicy {
  std::size_t base_intervals = 4096, maximum_fine_intervals = 65536;
  std::size_t maximum_total_work = 4000000;
  std::size_t maximum_native_bytes = 32 * 1024 * 1024;
  double absolute_fraction_tolerance = 1e-8, relative_fraction_tolerance = 2e-6;
  double absolute_opacity_tolerance = 2e-9, relative_opacity_tolerance = 3e-6;
};
struct HeliumEscapeHistoryWork {
  std::size_t imported_rows = 0, driver_evaluations = 0,
              rate_evaluations = 0;
  std::size_t total() const noexcept {
    return imported_rows + driver_evaluations + rate_evaluations;
  }
};
// Scope witnesses are separate from numerical error estimates and do not
// certify missing HeIII kinetics, atomic facts, collisions or the supplied law.
struct HeliumEscapeHistoryWitness {
  bool source_capture_complete = false;
  std::size_t attempted_rates = 0, nonfinite_temperature_attempts = 0;
  std::optional<double> minimum_attempted_temperature_kelvin,
      maximum_attempted_temperature_kelvin, minimum_attempted_tau,
      maximum_attempted_tau, maximum_log_retained_heiii_activity;
};
struct HeliumEscapeHistoryValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  // Empirical root/refinement/interpolation/arithmetic/cast diagnostic, not
  // a rigorous bound or physical input/model error.
  double absolute_error_estimate = 0;
};
inline constexpr unsigned helium_escape_fraction = 1,
                          helium_escape_electron_density = 2,
                          helium_escape_thomson_opacity = 4;
struct HeliumEscapeHistoryRow {
  double redshift = 0;
  HeliumEscapeHistoryValue helium_singly_ionized_fraction,
      electron_number_density_per_cubic_metre, thomson_opacity_per_redshift;
};
struct HeliumEscapeHistoryBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested_outputs = 0;
  std::size_t driver_evaluations = 0;
  std::vector<HeliumEscapeHistoryRow> rows;
};
class HeliumEscapeHistory {
public:
  HeliumEscapeHistory() = default;
  // Copies preflight simultaneous retained owners against the source policy;
  // resource refusal throws length_error, allocation failure throws bad_alloc.
  HeliumEscapeHistory(const HeliumEscapeHistory &);
  HeliumEscapeHistory &operator=(const HeliumEscapeHistory &);
  HeliumEscapeHistory(HeliumEscapeHistory &&) noexcept;
  HeliumEscapeHistory &operator=(HeliumEscapeHistory &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const HeliumEscapeDriver *source() const noexcept {
    return source_ ? &*source_ : nullptr;
  }
  HeliumEscapeHistoryWork work() const noexcept { return work_; }
  const HeliumEscapeHistoryWitness &witness() const noexcept { return witness_; }
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  // Queries preserve source order and individual refusals. No ODE/background
  // work; bounded geometric driver interpolation is charged to this batch.
  HeliumEscapeHistoryBatch evaluate(
      std::span<const double>, unsigned requested_outputs,
      std::size_t maximum_output_points = 4096,
      std::size_t maximum_native_bytes = 32 * 1024 * 1024) const;
private:
  struct Node {
    long double ell = 0, redshift = 0, fraction = 0,
                fraction_error = 0, curvature = 0;
  };
  numerics::Status status_ = numerics::Status::invalid_input;
  std::optional<HeliumEscapeDriver> source_;
  HeliumEscapeHistoryPolicy policy_;
  HeliumEscapeHistoryWork work_;
  HeliumEscapeHistoryWitness witness_;
  std::vector<Node> nodes_;
  friend HeliumEscapeHistory prepare_helium_escape_history(
      const HeliumEscapeHistoryRequest &, HeliumEscapeHistoryPolicy);
  friend std::optional<std::size_t> helium_escape_history_payload_bound(
      std::size_t, std::size_t, std::size_t, std::size_t) noexcept;
};
HeliumEscapeHistory prepare_helium_escape_history(
    const HeliumEscapeHistoryRequest &, HeliumEscapeHistoryPolicy = {});
// Conservative requested simultaneous preparation/output payload, with checked
// products. Actual vector/string capacities are checked after acquisition.
// Excludes borrowed inputs, allocator bookkeeping and RSS.
std::optional<std::size_t> helium_escape_history_payload_bound(
    std::size_t fine_intervals, std::size_t source_knots,
    std::size_t producer_identity_bytes, std::size_t output_points) noexcept;
} // namespace irred::cosmology
