// Conditional supplied linear power; numerical diagnostics are estimates, not bounds.
#pragma once
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace irred::windowed_linear_power {
inline constexpr char model_id[] =
    "supplied-linear-Kaiser-zero-dispersion-geometric-AP-window02/v1";
enum class Status {
  ok, invalid_input, nonfinite_input, outside_support, incompatible_epoch,
  unsupported_arithmetic, unresolved_projection, overflow, work_limit,
  payload_limit, numerical_budget_exceeded, invalid_owner
};
enum class Coordinates { physical_mpc, fixed_h_reference };
struct SourceIdentity {
  std::string spectrum, window, matter_component, growth_response;
  std::string grid_order, window_calibration;
};
struct Spectrum {
  double redshift;
  std::vector<double> k_per_mpc, power_mpc3;
};
struct Window {
  Coordinates coordinates;
  std::optional<double> h_reference; // engaged iff fixed_h_reference
  std::vector<double> theory_k, output0_k, output2_k;
  std::vector<std::uint64_t> output0_ids, output2_ids;
  std::vector<double> w00, w02, w20, w22; // explicit row-major names
};
struct Source { SourceIdentity identity; Spectrum spectrum; Window window; };
struct PreparationPolicy {
  std::size_t maximum_source_actions;
  std::size_t maximum_live_payload_bytes;
};
struct PreparationWork {
  std::size_t fields = 0, spectrum_nodes = 0, axis_nodes = 0;
  std::size_t window_coefficients = 0, mapped_axis_writes = 0;
};
struct ModelPoint {
  double redshift, b1, f, alpha_perp, alpha_parallel;
};
struct EvaluationPolicy {
  double maximum_multipole_absolute_error_estimate_mpc3;
  double maximum_mean_absolute_error_estimate_output_unit;
  std::size_t maximum_actions, maximum_angular_callbacks;
  std::size_t maximum_callbacks_per_integral;
  unsigned maximum_depth;
  std::size_t maximum_live_payload_bytes;
  bool retain_input_multipoles;
};
struct EvaluationWork {
  std::size_t model_fields = 0, support_columns = 0, split_candidates = 0;
  std::size_t angular_callbacks = 0, bracket_comparisons = 0;
  std::size_t interpolation_calls = 0, window_products = 0;
  std::size_t output_writes = 0;
};
struct PowerValue {
  double value = 0, quadrature_estimate = 0, arithmetic_estimate = 0;
  double grid_projection_estimate = 0, volume_projection_estimate = 0;
  double combined_estimate = 0;
};
struct MeanRow {
  unsigned multipole = 0;
  std::uint64_t source_row_id = 0;
  double source_k = 0;
  PowerValue power;
};
struct MultipoleRow {
  double source_k = 0;
  PowerValue p0, p2; // in the SAME output power unit as mean rows
};
struct Result {
  Status status = Status::invalid_owner;
  ModelPoint model{};
  SourceIdentity identity; // bounded metadata copy owns result provenance
  Coordinates coordinates = Coordinates::physical_mpc;
  std::optional<double> h_reference;
  PreparationWork preparation_work;
  EvaluationWork work; // all attempts, including refused callback actions
  std::size_t known_live_payload_estimate_bytes = 0;
  std::vector<MeanRow> mean; // all output0 then all output2, or empty on refusal
  std::vector<MultipoleRow> input_multipoles; // optional diagnostic output
  // Model adequacy/calibration/upstream uncertainty is unassessed, not zero.
  bool observational_qualified = false;
};
class Prepared {
public:
  Prepared() = default;
  Prepared(const Prepared&) = delete;
  Prepared& operator=(const Prepared&) = delete;
  Prepared(Prepared&&) noexcept;
  Prepared& operator=(Prepared&&) noexcept; // source invalidated; self no-op
  Status status() const noexcept;
  bool owns_source() const noexcept;
  const PreparationWork& preparation_work() const noexcept;
  std::optional<std::size_t> retained_payload_bytes() const noexcept;
  Result evaluate(const ModelPoint&, const EvaluationPolicy&) const;
private:
  Source source_{};
  std::vector<double> mapped_theory_k_per_mpc_;
  std::vector<long double> grid_projection_estimate_per_mpc_;
  long double maximum_abs_spectrum_slope_ = 0;
  Status status_ = Status::invalid_owner;
  PreparationWork work_;
  bool owns_source_ = false;
  friend Prepared prepare(Source&&, const PreparationPolicy&);
};
Prepared prepare(Source&&, const PreparationPolicy&);
} // namespace irred::windowed_linear_power
