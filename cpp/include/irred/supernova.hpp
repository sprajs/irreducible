#pragma once
#include "irred/background.hpp"
#include "irred/observations.hpp"
#include "irred/statistics.hpp"
#include <memory>
#include <span>
#include <variant>
#include <vector>
namespace irred::supernova {
enum class Status {
  ok,
  invalid_input,
  incompatible_metadata,
  numerical_failure,
  work_limit
};
struct MagnitudeCoordinate {
  double z_expansion;
  cosmology::Observer observer;
};
// Shared immutable source plus an explicitly declared selected/order-preserving
// view. Dataset cuts and release profile admission belong to named presets.
struct SelectedMagnitudeSource {
  std::shared_ptr<const observations::Prepared> source;
  std::vector<std::size_t> source_indices;
  std::vector<std::string> ordered_ids;
  std::vector<MagnitudeCoordinate> coordinates;
};
struct PreparationPolicy {
  numerics::Arithmetic arithmetic =
      static_cast<numerics::Arithmetic>(UINT32_MAX);
  std::size_t maximum_selected_rows = 0, maximum_matrix_elements = 0;
  double maximum_forward_sensitivity = 0;
  std::size_t maximum_native_bytes = 0;
};
// Conservative preparation peak, excluding preexisting shared source and
// allocator/RSS. Includes deterministic selected descriptor and metadata
// copies.
std::optional<std::size_t>
preparation_payload_bound(const observations::Prepared &, std::size_t,
                          numerics::Arithmetic) noexcept;
struct NoMagnitudeEffect {};
struct GreyLog1pMagnitude {
  double epsilon_mag;
  explicit GreyLog1pMagnitude(double e) : epsilon_mag(e) {}
};
using SourceEffectSpec = std::variant<NoMagnitudeEffect, GreyLog1pMagnitude>;
struct ModelPoint {
  cosmology::ExpansionSpec expansion;
  cosmology::FlatFLRW geometry;
  SourceEffectSpec source_effect;
  ModelPoint(cosmology::ExpansionSpec e, cosmology::FlatFLRW g,
             SourceEffectSpec s)
      : expansion(std::move(e)), geometry(g), source_effect(std::move(s)) {}
};
enum class Output : std::uint32_t {
  score = 1,
  geometric_shape = 2,
  magnitude_effect = 4,
  corrected_residuals = 8,
  profiled_residuals = 16,
  diagnostics = 32
};
struct Policy {
  cosmology::EvaluationPolicy background;
  numerics::Arithmetic arithmetic =
      static_cast<numerics::Arithmetic>(UINT32_MAX);
  std::size_t maximum_models = 0, maximum_native_bytes = 0;
  double maximum_forward_sensitivity = 0;
  std::uint32_t requested = 0;
};
struct Score {
  double offset_coefficient = 0, quadratic = 0, relative_profile_score = 0;
};
struct OutputState {
  cosmology::Availability availability = cosmology::Availability::not_requested;
  Status status = Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
};
struct Slot {
  ModelPoint source;
  Status status = Status::invalid_input;
  cosmology::Status background_status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  statistics::DensityStatus profile_status =
      statistics::DensityStatus::invalid_input;
  bool has_predictions = false;
  OutputState geometry, effect, corrected, profile;
  std::optional<Score> score;
  std::optional<statistics::ProfileResult> diagnostics;
  std::vector<double> geometric_shape, magnitude_shifts, corrected_residuals,
      profiled_residuals;
  std::vector<std::size_t> background_node_indices;
  cosmology::Work work;
  std::string_view model_id, hypothesis_id;
  explicit Slot(ModelPoint p) : source(std::move(p)) {}
};
struct BatchResult {
  Status status = Status::invalid_input;
  std::vector<Slot> slots;
  cosmology::Work work;
};
class Consumer {
public:
  Consumer() = default;
  Consumer(const Consumer &) = default;
  Consumer &operator=(const Consumer &) = default;
  // Distinct moves invalidate the source; self assignment retains its state.
  Consumer(Consumer &&) noexcept;
  Consumer &operator=(Consumer &&) noexcept;

  // Shared Prepared source excluded; context charges it once by owner identity.
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  Status status() const noexcept { return status_; }
  statistics::DensityStatus preparation_status() const noexcept {
    return preparation_status_;
  }
  numerics::Status preparation_numerical_status() const noexcept {
    return preparation_numerical_status_;
  }
  const SelectedMagnitudeSource &selected_source() const noexcept {
    return selected_;
  }
  const statistics::Metadata &probability_metadata() const noexcept {
    return profile_.metadata();
  }
  const statistics::ProfileOperator &profile_operator() const noexcept {
    return profile_;
  }
  static constexpr const char *score_id =
      "F03/unconstrained-offset-relative-profile/v1";
  static constexpr const char *prediction_unit = "magnitude";
  static constexpr const char *offset_convention =
      "unconstrained additive magnitude intercept; H0 absorbed";
  static constexpr const char *shape_convention =
      "declared expansion redshift radial integral with declared observer "
      "luminosity prefactor";
  BatchResult evaluate_batch(std::span<const ModelPoint>, Policy) const;

private:
  SelectedMagnitudeSource selected_;
  std::vector<double> observed_;
  statistics::ProfileOperator profile_;
  Status status_ = Status::invalid_input;
  statistics::DensityStatus preparation_status_ =
      statistics::DensityStatus::invalid_input;
  numerics::Status preparation_numerical_status_ =
      numerics::Status::invalid_input;
  friend Consumer prepare(SelectedMagnitudeSource, PreparationPolicy);
};
Consumer prepare(SelectedMagnitudeSource, PreparationPolicy);
// Named release adapter; no cut or table profile exists in the generic
// consumer.
SelectedMagnitudeSource
    pantheon_zhd_gt_001(std::shared_ptr<const observations::Prepared>);
} // namespace irred::supernova
