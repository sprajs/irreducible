#pragma once
#include "irred/growth_amplitude.hpp"
#include "irred/statistics.hpp"
namespace irred::rsd {
enum class RowRole : std::uint32_t { synthetic_control, unknown = UINT32_MAX };
enum class CovarianceUnit : std::uint32_t {
  dimensionless_f_sigma8_squared,
  unknown = UINT32_MAX
};
struct DensitySource {
  std::vector<double> scale_factors, observed;
  std::vector<std::string> ordered_ids, event_ids, covariance_axis_ids;
  RowRole role = RowRole::unknown;
  CovarianceUnit covariance_unit = CovarianceUnit::unknown;
  cosmology::Sigma8Convention amplitude_convention =
      cosmology::Sigma8Convention::unknown;
  std::string table_identity, covariance_identity, ordering_provenance,
      calibration_provenance, dependence_provenance;
};
struct DensityInput {
  DensitySource source;
  std::vector<double> covariance;
};
struct PreparationPolicy {
  size_t maximum_queries = 0, maximum_matrix_elements = 0,
         maximum_string_bytes = 0, maximum_native_bytes = 0;
  double maximum_forward_sensitivity = 0;
  numerics::Arithmetic arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
};
struct ModelPoint {
  cosmology::ExpansionSpec expansion;
  cosmology::FlatFLRW geometry;
  cosmology::FixedSigma8Source amplitude;
};
enum class Output : std::uint32_t {
  normalized_density = 1,
  predictions = 2,
  residuals = 4
};
struct OutputState {
  cosmology::Availability availability = cosmology::Availability::not_requested;
  numerics::Status numerical_status = numerics::Status::invalid_input;
};
struct DensityPolicy {
  cosmology::GrowthAmplitudePolicy predictions;
  size_t maximum_models = 0, maximum_queries = 0, maximum_string_bytes = 0,
         maximum_native_bytes = 0, maximum_total_callbacks = 0;
  double maximum_forward_sensitivity = 0;
  double maximum_projection_log_density_error = 1e-8;
  numerics::Arithmetic arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  unsigned requested = 0;
};
struct DensitySlot {
  ModelPoint source;
  explicit DensitySlot(ModelPoint x) : source(std::move(x)) {}
  OutputState predictions_state, residuals_state, density_state;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::optional<statistics::GaussianResult> result;
  std::vector<cosmology::GrowthAmplitudeValue> predictions, residuals;
  std::optional<cosmology::GrowthRow> reference;
  std::optional<double> projection_log_density_error_estimate;
  size_t callbacks = 0;
};
struct DensityBatch {
  statistics::DensityStatus status = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<DensitySlot> slots;
  size_t callbacks = 0;
};
std::optional<size_t>
retained_source_payload_bound(const DensitySource &) noexcept;
std::optional<size_t> preparation_payload_bound(const DensityInput &,
                                                numerics::Arithmetic) noexcept;
// Returned model sources/arrays plus serial growth, residual and solve scratch.
// strings includes max(32,size+1) for both copied amplitude strings per model.
std::optional<size_t> density_payload_bound(size_t models, size_t rows,
                                            size_t strings,
                                            unsigned requested) noexcept;
class PreparedDensity {
public:
  PreparedDensity() = default;
  PreparedDensity(const PreparedDensity &) = delete;
  PreparedDensity &operator=(const PreparedDensity &) = delete;
  PreparedDensity(PreparedDensity &&) noexcept;
  PreparedDensity &operator=(PreparedDensity &&) noexcept;
  statistics::DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  const DensitySource &source() const noexcept { return source_; }
  std::span<const double> covariance() const noexcept {
    return gaussian_.covariance();
  }
  const statistics::Metadata &metadata() const noexcept {
    return gaussian_.metadata();
  }
  std::optional<size_t> retained_payload_bound() const noexcept;
  DensityBatch evaluate(std::span<const ModelPoint>, DensityPolicy) const;

private:
  DensitySource source_;
  statistics::Gaussian gaussian_;
  statistics::DensityStatus status_ = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
  friend PreparedDensity prepare_density(DensityInput, PreparationPolicy);
};
PreparedDensity prepare_density(DensityInput, PreparationPolicy);
inline constexpr std::string_view conditional_density_id =
    "RSD/flat-pressureless-supplied-linear-fsigma8-synthetic-density/v1";
} // namespace irred::rsd
