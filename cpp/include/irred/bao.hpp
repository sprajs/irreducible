#pragma once
#include "irred/background.hpp"
#include "irred/statistics.hpp"
#include <optional>
#include <span>
#include <string>
#include <vector>
namespace irred::cosmology {
struct SoundHorizonRequest;
struct ThermalObservableRequest;
} // namespace irred::cosmology
namespace irred::bao {
struct ConditionalDensityPolicy;
struct ConditionalDensityBatch;
struct ThermalDensityPolicy;
struct ThermalDensityBatch;
enum class Observable : std::uint32_t {
  transverse_over_ruler,
  hubble_over_ruler,
  volume_over_ruler
};
struct Query {
  double z = 0;
  Observable observable = Observable::transverse_over_ruler;
};
struct Ruler {
  double h0_rd_km_s;
  explicit Ruler(double v) : h0_rd_km_s(v) {}
};
struct Policy {
  cosmology::EvaluationPolicy background;
  size_t maximum_queries = 0, maximum_native_bytes = 0;
};
struct Slot {
  Query source;
  cosmology::Status status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::optional<double> value;
  std::optional<size_t> node_index;
};
struct Batch {
  cosmology::Status status = cosmology::Status::invalid_input;
  std::vector<Slot> slots;
  cosmology::Work work;
};
inline constexpr std::string_view equation_id = "P01/flat-free-ruler-BAO/v1";
inline constexpr std::string_view ruler_convention_id =
    "P01/free-H0rd-km-s-no-early-physics/v1";
// E/I only; no computational H0 or observer coordinate. All covariance rows
// remain distinct while exact expansion-z nodes are reused by the provider.
Batch evaluate(const cosmology::Expansion &, Ruler, std::span<const Query>,
               Policy);
enum class RowRole : std::uint32_t {
  released_fitted_distance_summary = 0,
  synthetic_control = 1,
  unknown = UINT32_MAX
};
enum class CovarianceUnit : std::uint32_t {
  dimensionless_ratio_squared = 0,
  unknown = UINT32_MAX
};
struct DensityInput {
  std::vector<Query> queries;
  std::vector<double> observed, covariance;
  std::vector<std::string> ordered_ids;
  RowRole role = RowRole::unknown;
  CovarianceUnit covariance_unit = CovarianceUnit::unknown;
  std::string table_identity, covariance_identity, ordering_provenance,
      calibration_provenance, dependence_provenance;
};
// Transferred dynamic capacities, including conservative SSO string storage.
// Excludes object headers, allocator bookkeeping and RSS; overflow is absent.
std::optional<size_t>
retained_source_payload_bound(const DensityInput &) noexcept;
// Simultaneous preparation payload, including transferred source capacities
// and copied ID/provenance storage. Excludes borrowed input, allocator/RSS.
std::optional<size_t>
preparation_payload_bound(size_t rows, size_t source_payload_bytes,
                          size_t copied_identity_payload_bytes) noexcept;
struct PreparationPolicy {
  size_t maximum_queries = 0, maximum_matrix_elements = 0,
         maximum_string_bytes = 0, maximum_native_bytes = 0;
  double maximum_forward_sensitivity = 0;
  numerics::Arithmetic arithmetic =
      static_cast<numerics::Arithmetic>(UINT32_MAX);
};
struct ModelPoint {
  cosmology::ExpansionSpec expansion;
  cosmology::FlatFLRW geometry;
  Ruler ruler;
  ModelPoint(cosmology::ExpansionSpec e, cosmology::FlatFLRW g, Ruler r)
      : expansion(std::move(e)), geometry(g), ruler(r) {}
};
enum class Output : std::uint32_t {
  normalized_density = 1,
  predictions = 2,
  residuals = 4
};
struct DensityPolicy {
  Policy observables;
  size_t maximum_models = 0, maximum_native_bytes = 0;
  double maximum_forward_sensitivity = 0;
  numerics::Arithmetic arithmetic =
      static_cast<numerics::Arithmetic>(UINT32_MAX);
  std::uint32_t requested = 0;
};
struct OutputState {
  cosmology::Availability availability = cosmology::Availability::not_requested;
  cosmology::Status status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
};
struct DensitySlot {
  ModelPoint source;
  explicit DensitySlot(ModelPoint p) : source(std::move(p)) {}
  std::string_view model_id;
  cosmology::Status background_status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  bool has_predictions = false;
  OutputState predictions_state, residuals_state, density_state;
  std::optional<statistics::GaussianResult> result;
  std::vector<double> predictions, residuals;
  std::vector<size_t> node_indices;
  cosmology::Work work;
};
struct DensityBatch {
  statistics::DensityStatus status = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<DensitySlot> slots;
  cosmology::Work work;
};
class PreparedDensity {
public:
  PreparedDensity() = default;
  PreparedDensity(const PreparedDensity &) = default;
  PreparedDensity &operator=(const PreparedDensity &) = default;
  // Distinct moves invalidate the source; self assignment retains its state.
  PreparedDensity(PreparedDensity &&) noexcept;
  PreparedDensity &operator=(PreparedDensity &&) noexcept;

  std::optional<std::size_t> retained_payload_bound() const noexcept;
  statistics::DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  const DensityInput &source() const noexcept { return source_; }
  const statistics::Metadata &metadata() const noexcept {
    return gaussian_.metadata();
  }
  DensityBatch evaluate(std::span<const ModelPoint>, DensityPolicy) const;
  ConditionalDensityBatch
      evaluate_conditional(std::span<const cosmology::SoundHorizonRequest>,
                           ConditionalDensityPolicy) const;

  ThermalDensityBatch
      evaluate_thermal(std::span<const cosmology::ThermalObservableRequest>,
                       ThermalDensityPolicy) const;

private:
  DensityInput source_;
  statistics::Gaussian gaussian_;
  statistics::DensityStatus status_ = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
  friend PreparedDensity prepare_density(DensityInput, PreparationPolicy);
};
// Explicit owning input transfer; same full normalized covariance measure.
PreparedDensity prepare_density(DensityInput, PreparationPolicy);
} // namespace irred::bao
